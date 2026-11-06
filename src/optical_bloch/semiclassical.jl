# Semiclassical cooling rates: the internal master equation under the Doppler shifts
# of classical coherent-state motion, with the rate of energy change from the
# time-averaged laser force (the "PACMAN" method of [Bartolotta2024]).

"""
A numerical engine of the semiclassical cooling rates:
[`HarmonicBalance`](@ref) or [`TimeAverage`](@ref).
"""
abstract type SemiclassicalEngine end

"""
    HarmonicBalance(; margin = 16)

Solves the periodic steady state of a single oscillating mode directly, as a
block-tridiagonal system over the harmonics ``ρ_q`` of the mode frequency: there
is no transient and no averaging window. Its cost grows linearly with the number
of harmonics kept, ``⌈β⌉ +`` `margin`, with ``β`` the largest difference of Doppler
modulation indices between two states. The cutoff is raised automatically
if the last harmonic is not negligible.
"""
struct HarmonicBalance <: SemiclassicalEngine
    margin::Int
end

function HarmonicBalance(; margin::Integer=16)
    margin >= 1 || throw(ArgumentError("margin must be positive, got $margin"))
    HarmonicBalance(margin)
end

"""
    TimeAverage(; settling_time = 5u"µs", averaging_time = 20u"µs", num_phases = 16, window = :hann, reltol = 1e-6, abstol = 1e-8)

Integrates the master equation along the classical trajectory and averages the
power over `averaging_time` after `settling_time`, starting from the steady
state without motion, and over `num_phases` random sets of secular phases. This
is the only engine for several oscillating modes, whose incommensurate frequencies
have no common period. Both times and `num_phases` need convergence checks: the
averaging time must far exceed the inverse of the smallest mode frequency
difference, which the default only suits for differences above a few 100 kHz (angular
frequency difference × time ≫ 2π), or the beat of two modes adds to their rates. `window = :hann` weights the average by a Hann window, suppressing the
beats between modes much better than the `:rectangular` average of the original
method at equal length.
"""
struct TimeAverage{T<:Quantity} <: SemiclassicalEngine
    settling_time::T
    averaging_time::T
    num_phases::Int
    window::Symbol
    reltol::Float64
    abstol::Float64
end

function TimeAverage(;
    settling_time=5.0u"µs",
    averaging_time=20.0u"µs",
    num_phases::Integer=16,
    window::Symbol=:hann,
    reltol::Real=1e-6,
    abstol::Real=1e-8,
)
    settling_time >= zero(settling_time) ||
        throw(ArgumentError("settling_time must not be negative"))
    averaging_time > zero(averaging_time) ||
        throw(ArgumentError("averaging_time must be positive"))
    num_phases >= 1 || throw(ArgumentError("num_phases must be positive"))
    window in (:hann, :rectangular) ||
        throw(ArgumentError("window must be :hann or :rectangular, got $window"))
    settling, averaging = promote(float(settling_time), float(averaging_time))
    TimeAverage(
        settling,
        averaging,
        Int(num_phases),
        window,
        Float64(reltol),
        Float64(abstol),
    )
end

"""
    SemiclassicalRates(; occupation, engine = nothing)

Selects the semiclassical rates of [`semiclassical_cooling_rates`](@ref) for
[`cooling_metrics`](@ref), at the given coherent-state `occupation` (one value or
one per mode). `τ_c` is the local ``1/W_\\mathrm{SC}`` at that occupation (negative
if the mode is heated), not the area of a transient. `nbar` is `NaN`: without recoil,
the model has no steady occupation. Accepts no `heating_rate`.
"""
struct SemiclassicalRates{E} <: CoolingMethod
    occupation::Union{Float64,Vector{Float64}}
    engine::E
end

function SemiclassicalRates(;
    occupation,
    engine::Union{Nothing,SemiclassicalEngine}=nothing,
)
    all(>(0), occupation) || throw(ArgumentError("occupation must be positive"))
    SemiclassicalRates(
        occupation isa Real ? Float64(occupation) : collect(Float64, occupation),
        engine,
    )
end

# The plain-number form of a model and its modes, on column-stacked density matrices.
struct SemiclassicalProblem{T<:Unitful.Units}
    num_states::Int
    "Static Liouvillian, in the frame without motion."
    liouvillian::SparseMatrixCSC{ComplexF64,Int}
    "Positions of the populations in the column-stacked density matrix."
    diagonal::Vector{Int}
    "Doppler operators K_μ (the diagonal, state × mode)."
    doppler::Matrix{Float64}
    "Mode frequencies."
    frequencies::Vector{Float64}
    "Lamb–Dicke factors (beam × mode)."
    lamb_dicke_factors::Matrix{Float64}
    "Per beam, the indices and values of the row a with a ⋅ vec(ρ) = Tr[i(C - C†)/2 ρ]."
    force_indices::Vector{Vector{Int}}
    force_values::Vector{Vector{ComplexF64}}
    time_unit::T
end

# The Doppler gauge removes the beam phases only if each coupling crosses exactly its
# own beam once, and leaves the decay operators unchanged up to a phase only if all
# their components see one phase; both hold for a static frame by construction but for
# accidental exact beat-free loops.
function check_doppler_gauge(model::LindbladModel)
    traversals = model.frame.traversals
    num_beams = size(traversals, 2)
    for (b, C) in enumerate(model.couplings)
        scale = maximum(abs, C; init=zero(real(eltype(C))))
        for index in CartesianIndices(C)
            abs(C[index]) > 1e-12 * scale || continue
            upper, lower = Tuple(index)
            step = traversals[upper, :] .- traversals[lower, :]
            all(step[j] == (j == b) for j in 1:num_beams) || throw(
                ArgumentError(
                    "Beam $b closes a loop of beams on the rotating-frame tree, so its " *
                    "Doppler phase cannot be moved into the state energies",
                ),
            )
        end
    end
    for L in first(decay_operators(model))
        steps = [
            traversals[Tuple(index)[1], :] .- traversals[Tuple(index)[2], :] for
            index in CartesianIndices(L) if !iszero(L[index])
        ]
        all(==(first(steps)), steps) || throw(
            ArgumentError(
                "A decay operator connects states with different beam histories; " *
                "use decay = :resolved",
            ),
        )
    end
end

function SemiclassicalProblem(
    model::LindbladModel,
    mcs::AbstractVector{<:MotionalCoupling};
    time_unit::Unitful.Units=u"µs",
)
    require_static(model, "The semiclassical rate")
    isempty(mcs) && throw(ArgumentError("At least one mode is required"))
    check_doppler_gauge(model)
    num_states = length(model.basis)
    H = sparse(ComplexF64.(ustrip.(time_unit^-1, model.hamiltonian)))
    jumps = [
        sparse(ComplexF64.(ustrip.(time_unit^(-1 // 2), L))) for
        L in model.jump_operators
    ]
    L = -im * commutator_superoperator(H)
    for c in jumps
        L += dissipator_superoperator(c)
    end
    factors = reduce(hcat, [mc.projected_lamb_dicke for mc in mcs])
    force_rows = map(model.couplings) do C
        C_stripped = ComplexF64.(ustrip.(time_unit^-1, C))
        # Tr[G ρ] = ∑_rc G_cr ρ_rc
        row = vec(transpose(im .* (C_stripped .- C_stripped') ./ 2))
        scale = maximum(abs, row; init=0.0)
        indices = findall(x -> abs(x) > 1e-14 * scale, row)
        indices, row[indices]
    end
    SemiclassicalProblem(
        num_states,
        L,
        collect(1:(num_states+1):(num_states^2)),
        Float64.(model.frame.traversals) * factors,
        [ustrip(time_unit^-1, mc.mode.frequency) for mc in mcs],
        factors,
        first.(force_rows),
        last.(force_rows),
        time_unit,
    )
end

# The value a ⋅ u of each beam's force row.
function beam_forces!(forces, p::SemiclassicalProblem, u)
    for b in eachindex(forces)
        indices, values = p.force_indices[b], p.force_values[b]
        total = 0.0
        @inbounds for j in eachindex(indices)
            total += real(values[j] * u[indices[j]])
        end
        forces[b] = total
    end
    forces
end

function force_row(p::SemiclassicalProblem, mode)
    row = zeros(ComplexF64, p.num_states^2)
    for b in eachindex(p.force_indices)
        row[p.force_indices[b]] .+= p.lamb_dicke_factors[b, mode] .* p.force_values[b]
    end
    row
end

# Per column-stacked density-matrix element, K_r - K_c of one mode.
function doppler_differences(p::SemiclassicalProblem, mode)
    K = view(p.doppler, :, mode)
    vec([K[r] - K[c] for r in 1:p.num_states, c in 1:p.num_states])
end

# --- One oscillating mode: the harmonics of its periodic steady state ----------------

# Returns ⟨ṅ⟩ of the mode at occupation n. With ρ(t) = ∑_q ρ_q e^{-iqωt}, the equations
# (L + iqω) ρ_q + (A/2) d (ρ_{q+1} - ρ_{q-1}) = 0 close downwards from the cutoff on
# ρ_q = S_q ρ_{q-1}, and ρ_{-q} = ρ_q†.
function harmonic_balance_rate(
    p::SemiclassicalProblem,
    mode::Integer,
    n::Real,
    engine::HarmonicBalance,
)
    ω = p.frequencies[mode]
    d = doppler_differences(p, mode)
    # ρ†: a transposition of the density-matrix indices with conjugation.
    N = p.num_states
    adjoint_order = vec(permutedims(reshape(1:(N^2), N, N)))
    half_drive = sqrt(n) * ω .* d
    modulation = maximum(abs, 2 * sqrt(n) .* d; init=0.0)
    L = Matrix(p.liouvillian)
    margin = engine.margin
    for attempt in 1:4
        cutoff = ceil(Int, modulation) + margin
        ρ0, ρ1, tail = harmonic_recursion(p, L, half_drive, ω, cutoff, adjoint_order)
        if tail <= 1e-9 * norm(ρ0) || attempt == 4
            tail <= 1e-6 * norm(ρ0) || @warn(
                "Harmonic balance not converged: last harmonic at $(tail / norm(ρ0)) " *
                "of the mean for n = $n"
            )
            # R = 2√n Im g₁: the sin(ωt) component of g(t) = ∑_q g_q e^{-iqωt}.
            return 2 * sqrt(n) * imag(dot(conj(force_row(p, mode)), ρ1))
        end
        margin *= 2
    end
end

function harmonic_recursion(p, L, half_drive, ω, cutoff, adjoint_order)
    num_entries = size(L, 1)
    D = Diagonal(half_drive)
    steps = Vector{Matrix{ComplexF64}}(undef, cutoff)
    M = Matrix{ComplexF64}(undef, num_entries, num_entries)
    for q in cutoff:-1:1
        copyto!(M, L)
        q < cutoff && (M .+= D * steps[q+1])
        @inbounds for i in 1:num_entries
            M[i, i] += im * q * ω
        end
        steps[q] = lu!(M) \ Matrix(D)
    end
    S = steps[1]
    # ρ_{-1} = (S ρ₀)† = T ρ₀, with ρ₀ Hermitian.
    T = conj.(S[adjoint_order, adjoint_order])
    G = sparse(L .+ D * (S .- T))
    ρ0 = stationary_vector(G, p.diagonal)
    ρ = ρ0
    ρ1 = steps[1] * ρ0
    for step in steps
        ρ = step * ρ
    end
    ρ0, ρ1, norm(ρ)
end

# --- Several oscillating modes: time averaging along the trajectory -------------------

mutable struct TrajectoryState{P<:SemiclassicalProblem}
    problem::P
    amplitudes::Vector{Float64}
    phases::Vector{Float64}
    forces::Vector{Float64}
    state_shifts::Vector{Float64}
    drives::Vector{Float64}
    window_start::Float64
    window_length::Float64
    hann::Bool
end

# In-place ρ̇ = L ρ - i [∑_μ -A_μ sin(ω_μ t + θ_μ) K_μ, ρ].
function trajectory_liouvillian!(du, u, s::TrajectoryState, t)
    p = s.problem
    N = p.num_states
    fill!(du, zero(eltype(du)))
    L = p.liouvillian
    rows, values = rowvals(L), nonzeros(L)
    @inbounds for column in 1:(N^2)
        x = u[column]
        for j in nzrange(L, column)
            du[rows[j]] += values[j] * x
        end
    end
    @inbounds for μ in eachindex(s.drives)
        s.drives[μ] = sin(p.frequencies[μ] * t + s.phases[μ])
    end
    fill!(s.state_shifts, 0.0)
    @inbounds for μ in eachindex(s.drives), i in 1:N
        s.state_shifts[i] -=
            s.amplitudes[μ] * p.frequencies[μ] * s.drives[μ] * p.doppler[i, μ]
    end
    k = 0
    @inbounds for c in 1:N, r in 1:N
        k += 1
        du[k] -= im * (s.state_shifts[r] - s.state_shifts[c]) * u[k]
    end
    nothing
end

# The right-hand side with the power accumulators after the density matrix.
function averaging_rhs!(du, u, s::TrajectoryState, t)
    p = s.problem
    trajectory_liouvillian!(du, u, s, t)
    beam_forces!(s.forces, p, u)
    x = (t - s.window_start) / s.window_length
    weight = s.hann ? 2 * sinpi(x)^2 : 1.0
    offset = p.num_states^2
    @inbounds for μ in eachindex(s.drives)
        force = 0.0
        for b in eachindex(s.forces)
            force += p.lamb_dicke_factors[b, μ] * s.forces[b]
        end
        du[offset+μ] = weight * s.amplitudes[μ] * s.drives[μ] * force
    end
    nothing
end

# Returns ⟨ṅ_μ⟩ per mode of one trajectory with amplitudes 2√n_μ.
function trajectory_rates(
    p::SemiclassicalProblem,
    occupations,
    phases,
    engine::TimeAverage,
)
    num_entries, num_modes = p.num_states^2, length(p.frequencies)
    settling = ustrip(p.time_unit, engine.settling_time)
    averaging = ustrip(p.time_unit, engine.averaging_time)
    s = TrajectoryState(
        p,
        2 .* sqrt.(occupations),
        collect(Float64, phases),
        zeros(length(p.force_indices)),
        zeros(p.num_states),
        zeros(num_modes),
        settling,
        averaging,
        engine.window == :hann,
    )
    options = (;
        reltol=engine.reltol,
        abstol=engine.abstol,
        save_everystep=false,
        save_start=false,
    )
    ρ = stationary_vector(p.liouvillian, p.diagonal)
    if settling > 0
        settled = solve(
            ODEProblem(
                (du, u, q, t) -> trajectory_liouvillian!(du, u, q, t),
                ρ,
                (0.0, settling),
                s,
            ),
            Vern7();
            options...,
        )
        check_trajectory(settled)
        ρ = last(settled.u)
    end
    u0 = [ρ; zeros(ComplexF64, num_modes)]
    averaged = solve(
        ODEProblem(
            (du, u, q, t) -> averaging_rhs!(du, u, q, t),
            u0,
            (settling, settling + averaging),
            s,
        ),
        Vern7();
        options...,
    )
    check_trajectory(averaged)
    real.(last(averaged.u)[(num_entries+1):end]) ./ averaging
end

function check_trajectory(solution)
    successful_retcode(solution) ||
        error("Semiclassical trajectory failed: $(solution.retcode)")
end

# --- Front-ends ------------------------------------------------------------------------

function check_occupations(occupations, num_modes)
    length(occupations) == num_modes ||
        throw(ArgumentError("occupations must be one value per mode"))
    all(n -> isfinite(n) && n >= 0, occupations) ||
        throw(ArgumentError("occupations must be non-negative, got $occupations"))
end

default_engine(occupations) =
    count(>(0), occupations) <= 1 ? HarmonicBalance() : TimeAverage()

# ⟨ṅ_μ⟩ per mode, for the occupations of one sample and one engine.
function phonon_rates(p::SemiclassicalProblem, occupations, engine::HarmonicBalance, _)
    hot = findall(>(0), occupations)
    length(hot) <= 1 || throw(
        ArgumentError(
            "HarmonicBalance handles one oscillating mode, got $(length(hot)); use TimeAverage",
        ),
    )
    rates = zeros(length(occupations))
    for μ in hot
        rates[μ] = harmonic_balance_rate(p, μ, occupations[μ], engine)
    end
    rates
end

function phonon_rates(p::SemiclassicalProblem, occupations, engine::TimeAverage, rng)
    phases = 2π .* rand(rng, length(occupations), engine.num_phases)
    results = Vector{Vector{Float64}}(undef, engine.num_phases)
    Threads.@threads for j in 1:(engine.num_phases)
        results[j] = trajectory_rates(p, occupations, view(phases, :, j), engine)
    end
    reduce(+, results) ./ engine.num_phases
end

"""
    semiclassical_phonon_rates(model::LindbladModel, mcs::AbstractVector{<:MotionalCoupling}, occupations; engine = nothing, rng = default_rng(), time_unit = u"µs")

Returns the mean rates of change ``R_μ = ⟨\\dot n_μ⟩`` (inverse `time_unit`) of the
occupations of the given coherent-state motion:
``⟨\\hat a_μ⟩ = \\sqrt{n_μ} e^{-i(ω_μ t + θ_μ)}``, held at fixed energy, with the
secular phases ``θ_μ`` averaged out. This is the semiclassical method of
[Bartolotta2024](@citet) (PACMAN). Replacing the motion by its expectation value turns
each beam's displacement operator into the phase
``φ_b(t) = ∑_μ 2η_{bμ}\\sqrt{n_μ} \\cos(ω_μ t + θ_μ)``, and a transformation with the
beams' traversals of the rotating-frame tree moves those phases into Doppler shifts,

```math
H(t) = H - ∑_μ 2\\sqrt{n_μ} ω_μ \\sin(ω_μ t + θ_μ) K_μ,
\\qquad K_μ = \\mathrm{diag}\\Big(∑_b n_{i,b} η_{bμ}\\Big),
\\qquad
R_μ = \\Big⟨2\\sqrt{n_μ} \\sin(ω_μ t + θ_μ) \\mathrm{Tr}[H_{\\mathrm{sb},μ} ρ(t)]\\Big⟩,
```

with `H` the model's static Hamiltonian and ``H_{\\mathrm{sb},μ}`` the mode's
`sideband_hamiltonian`.

The expansion in ``η`` and the separation of timescales are both avoided, and the
modes are coupled through their common Doppler signal, so each ``R_μ`` depends on the
energy in all modes. Conversely, there is no recoil and no heating from the blue
sideband (``A_+``), hence no steady occupation: the rates are the drift of a hot mode,
reliable for ``η\\sqrt{2n+1} ≳ 0.1``, and in the limit ``n → 0`` return ``-n (A_- - A_+)``
of [`AdiabaticElimination`](@ref).

`engine` defaults to [`HarmonicBalance`](@ref) for at most one non-zero occupation and
to [`TimeAverage`](@ref), with random phases from `rng`, otherwise. Requires a static
model. Rates of modes with zero occupation are zero.
"""
function semiclassical_phonon_rates(
    model::LindbladModel,
    mcs::AbstractVector{<:MotionalCoupling},
    occupations;
    engine::Union{Nothing,SemiclassicalEngine}=nothing,
    rng::AbstractRNG=default_rng(),
    time_unit::Unitful.Units=u"µs",
)
    check_occupations(occupations, length(mcs))
    p = SemiclassicalProblem(model, mcs; time_unit)
    phonon_rates(p, occupations, something(engine, default_engine(occupations)), rng) .*
    time_unit^-1
end

"""
    semiclassical_cooling_rates(model::LindbladModel, mcs::AbstractVector{<:MotionalCoupling}, occupations; engine = nothing, rng = default_rng(), time_unit = u"µs")

Returns the energy-dependent cooling rates ``W_μ = -R_μ / n_μ`` per mode (inverse
`time_unit`; negative where heated, `NaN` for zero occupation) of
[`semiclassical_phonon_rates`](@ref).
"""
semiclassical_cooling_rates(
    model::LindbladModel,
    mcs::AbstractVector{<:MotionalCoupling},
    occupations;
    kwargs...,
) = .-semiclassical_phonon_rates(model, mcs, occupations; kwargs...) ./ occupations

"""
    thermal_cooling_rates(model::LindbladModel, mcs::AbstractVector{<:MotionalCoupling}, nbars; num_samples = 100, engine = nothing, rng = default_rng(), time_unit = u"µs") -> (; rates, standard_errors)

Returns the cooling rates ``W_{\\mathrm{th},μ} = -⟨R_μ⟩/\\bar n_μ`` of a product of
thermal states of the mean occupations `nbars` per mode (inverse `time_unit`), and the
standard errors of this Monte Carlo estimate from `num_samples` draws of the coherent
occupations, whose distribution is the Glauber–Sudarshan ``P`` function: independent
exponentials of means `nbars` [Bartolotta2024; Eqs. (31)–(34)](@cite). The estimate is
the ratio of the sample means of ``R_μ`` and ``n_μ``, which removes the noise of the
drawn occupations where ``R_μ ∝ n_μ``, with a bias of order 1/`num_samples`. `engine`
defaults to [`HarmonicBalance`](@ref) for one mode and [`TimeAverage`](@ref) otherwise.
Every random number is drawn from `rng` before the threaded evaluation, so results do
not depend on the number of threads.
"""
function thermal_cooling_rates(
    model::LindbladModel,
    mcs::AbstractVector{<:MotionalCoupling},
    nbars;
    num_samples::Integer=100,
    engine::Union{Nothing,SemiclassicalEngine}=nothing,
    rng::AbstractRNG=default_rng(),
    time_unit::Unitful.Units=u"µs",
)
    num_samples >= 2 || throw(ArgumentError("num_samples must be at least 2"))
    num_modes = length(mcs)
    check_occupations(nbars, num_modes)
    all(>(0), nbars) || throw(ArgumentError("nbars must be positive, got $nbars"))
    p = SemiclassicalProblem(model, mcs; time_unit)
    selected = something(engine, num_modes == 1 ? HarmonicBalance() : TimeAverage())
    samples = randexp(rng, num_modes, num_samples) .* collect(float.(nbars))
    num_phases = selected isa TimeAverage ? selected.num_phases : 1
    phases = 2π .* rand(rng, num_modes, num_phases, num_samples)
    rates = zeros(num_modes, num_phases, num_samples)
    Threads.@threads for index in CartesianIndices((num_phases, num_samples))
        j, k = Tuple(index)
        rates[:, j, k] .=
            selected isa TimeAverage ?
            trajectory_rates(p, view(samples, :, k), view(phases, :, j, k), selected) :
            phonon_rates(p, view(samples, :, k), selected, rng)
    end
    per_sample = dropdims(sum(rates; dims=2); dims=2) ./ num_phases
    # The ratio of the sample means: R ∝ n at low energy, so that the noise of the drawn
    # occupations themselves cancels, as it would not on dividing by the exact mean.
    mean_occupations = vec(sum(samples; dims=2)) ./ num_samples
    W = .-vec(sum(per_sample; dims=2)) ./ num_samples ./ mean_occupations
    residuals = per_sample .+ W .* samples
    variances = vec(sum(abs2, residuals; dims=2)) ./ (num_samples - 1)
    (;
        rates=W .* time_unit^-1,
        standard_errors=(sqrt.(variances ./ num_samples) ./ mean_occupations) .*
                        time_unit^-1,
    )
end

function cooling_metrics(
    model::LindbladModel,
    mcs::AbstractVector{<:MotionalCoupling},
    method::SemiclassicalRates;
    heating_rate=nothing,
    ρ=nothing,
    time_unit::Unitful.Units=u"µs",
)
    isnothing(heating_rate) || throw(
        ArgumentError(
            "The semiclassical rates describe only the laser drift, not a heating rate",
        ),
    )
    num_modes = length(mcs)
    occupations =
        method.occupation isa Real ? fill(method.occupation, num_modes) :
        method.occupation
    rates = semiclassical_cooling_rates(
        model,
        mcs,
        occupations;
        engine=method.engine,
        time_unit,
    )
    [
        CoolingMetrics(method, NaN, 1 / W, nothing, nothing, nothing, nothing) for
        W in rates
    ]
end

export HarmonicBalance, TimeAverage, SemiclassicalRates
export semiclassical_phonon_rates, semiclassical_cooling_rates, thermal_cooling_rates
