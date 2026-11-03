# Nearest-neighbour phonon rate equations that stay valid for saturated sidebands:
# their rates are evaluated non-perturbatively in the sideband coupling, rung by
# rung. Rung n joins the Fock states n - 1 and n.

"""
    SaturatedRates(; initial_nbar, num_fock = fock_truncation.(initial_nbar), lamb_dicke_order = Inf, bandwidth = 2, node_ratio = 2)

Selects the steady occupation and integral relaxation time of a
[`PhononRateEquation`](@ref) for [`cooling_metrics`](@ref). `τ_c` has the same
meaning as for [`IntegratedTransient`](@ref), from a thermal state of
`initial_nbar` (one value or one per mode), for a small fraction of its cost.
`num_fock`, `lamb_dicke_order`, `bandwidth` and `node_ratio` pass to
[`PhononRateEquation`](@ref). Modes are treated independently. The result
carries the rate equation for cooling curves and other initial states.
"""
struct SaturatedRates <: CoolingMethod
    initial_nbar::Union{Float64,Vector{Float64}}
    num_fock::Union{Nothing,Int,Vector{Int}}
    lamb_dicke_order::Float64
    bandwidth::Int
    node_ratio::Float64
end

function SaturatedRates(;
    initial_nbar,
    num_fock=nothing,
    lamb_dicke_order::Real=Inf,
    bandwidth::Integer=2,
    node_ratio::Real=2,
)
    validate_rate_settings(lamb_dicke_order, bandwidth, node_ratio)
    all(>=(0), initial_nbar) ||
        throw(ArgumentError("initial_nbar must be non-negative"))
    SaturatedRates(
        initial_nbar isa Real ? Float64(initial_nbar) : collect(Float64, initial_nbar),
        isnothing(num_fock) ? nothing :
        num_fock isa Integer ? Int(num_fock) : collect(Int, num_fock),
        lamb_dicke_order,
        bandwidth,
        node_ratio,
    )
end

function validate_rate_settings(lamb_dicke_order, bandwidth, node_ratio)
    lamb_dicke_order == 1 ||
        isinf(lamb_dicke_order) ||
        throw(ArgumentError("lamb_dicke_order must be 1 or Inf, got $lamb_dicke_order"))
    bandwidth >= 1 || throw(ArgumentError("bandwidth must be positive, got $bandwidth"))
    node_ratio >= 1 || throw(ArgumentError("node_ratio must be at least 1"))
    nothing
end

"""
    PhononRateEquation(model::LindbladModel, mc::MotionalCoupling; num_fock, lamb_dicke_order = Inf, bandwidth = 2, node_ratio = 2, heating_rate = 0u"s^-1", time_unit = u"µs")

The nearest-neighbour rate equations ``\\dot p_n = Γ^↑_n p_{n-1} + Γ^↓_{n+1}
p_{n+1} - (Γ^↓_n + Γ^↑_{n+1}) p_n`` of the phonon populations of one mode, with
`lowering_rates` ``Γ^↓_n`` (``n → n - 1``) and `raising_rates` ``Γ^↑_n``
(``n - 1 → n``) for ``n = 1 … \\mathtt{num\\_fock} - 1``. Unlike
[`AdiabaticElimination`](@ref), whose rates ``nA_∓`` are second order in the
sideband coupling, these saturate once ``ηΩ\\sqrt{n}`` reaches the internal
relaxation rate.

Fock-resolved rate equations are the low-saturation result of Morigi et al.,
Phys. Rev. A **59**, 3797 (1999), Eq. (15), and saturated sideband transfer in
a two-state model that of Phatak et al., Phys. Rev. A **110**, 043116 (2024),
Eqs. (66)–(68). Here, each rung is solved for arbitrary internal structure and
drive, from the steady state of a ladder with the couplings of rung ``n`` on
every rung: a Floquet problem with `bandwidth` harmonics, equivalent to the
power absorbed from an oscillation of amplitude ``∝\\sqrt{n}``. The bare
population of a Fock state would not do: on a saturated sideband, it is
coherently shared with its neighbour. The net rate ``Γ^↑_n - Γ^↓_n`` is the
mean phonon-number current of the ladder, and the total ``Γ^↑_n + Γ^↓_n`` its
phonon-number diffusion: the second cumulant of its counting statistics, with
the recoil of the emissions. Where the steps are directional (``Γ^↑_n/Γ^↓_n ≲
0.8``, as they are where a cold steady occupation is decided) and the sideband
reach ``2|η|\\sqrt{n}`` stays below one, the ratio ``Γ^↑_n/Γ^↓_n`` follows
instead from the steady state of the closed internal ⊗ Fock system of the
states ``n-2 … n+1``, which gives it directly rather than as the small
difference of two rates; so it does wherever the ladder's steps are
sub-Poissonian (a diffusion below the net rate, on a saturated sideband),
which nearest-neighbour steps cannot be.

For weak sidebands, the rates are ``Γ^↓_n = nA_-`` and ``Γ^↑_n = nA_+``.
`lamb_dicke_order = Inf` uses the full
displacement elements of [`displacement_elements`](@ref) in the couplings (the
carrier averaged over the rung) and the emission rule of [`emission_rule`](@ref)
in the recoil; `1` follows the first-order [`motional_model`](@ref) and also
accepts a constant `recoil_moment`. Rungs are solved at nodes at most
`node_ratio` apart, refined around zeros of the carrier and first-sideband
elements (dark Fock states), and interpolated between (the net rate per
quantum linearly in ``\\ln n``, the slower rate per quantum in logs);
`node_ratio = 1` solves every rung. `heating_rate` adds ``n \\dot
n_\\mathrm{heat}`` to both rates.

As independent steps between neighbouring Fock states, this cannot represent
sub-Poissonian statistics or the bypassing of a dark Fock state through a second
sideband, next to which the ladder also understates the diffusion. Around a
reach of one, both sources of the ratio are off by tens of percent, which
affects the width of the distribution but neither its mean nor a cold steady
occupation. Modes are independent: their competition for the same scattering
events is not captured. Requires an internal model without beat harmonics. A
heated mode accumulates at the truncation, where `steady_value` depends on
`num_fock`.

[`steady_value`](@ref), [`fock_populations`](@ref),
[`integral_relaxation_time`](@ref), [`cooling_curve`](@ref) and
[`time_to_nbar`](@ref) evaluate the equation.
"""
struct PhononRateEquation{R<:Quantity,T<:Unitful.Units}
    lowering_rates::Vector{R}
    raising_rates::Vector{R}
    time_unit::T
end

function PhononRateEquation(
    model::LindbladModel,
    mc::MotionalCoupling;
    num_fock::Integer,
    lamb_dicke_order::Real=Inf,
    bandwidth::Integer=2,
    node_ratio::Real=2,
    heating_rate=0.0u"s^-1",
    time_unit::Unitful.Units=u"µs",
)
    require_static(model, "PhononRateEquation")
    validate_rate_settings(lamb_dicke_order, bandwidth, node_ratio)
    num_fock >= 2 || throw(ArgumentError("num_fock must be at least 2, got $num_fock"))
    heating = ustrip(time_unit^-1, only(mode_heating_rates(heating_rate, 1)))
    ladder = RungLadder(model, mc; lamb_dicke_order, bandwidth, time_unit)
    rungs = 1:(num_fock-1)
    nodes = rung_nodes(ladder, last(rungs), node_ratio)
    rates = [rung_rates(ladder, n) for n in nodes]
    lowering, raising = first.(rates), last.(rates)
    # The rates grow ∝ n without saturation, hence interpolation of the rates per
    # quantum. The net rate is the physical quantity and nearly cancels at high n,
    # so it is interpolated as such, and the slower rate in logs, as it spans decades
    # across a dark Fock state. This keeps both rates positive.
    if nodes != rungs
        net = linear_interpolation(log.(nodes), (raising .- lowering) ./ nodes)
        slower =
            linear_interpolation(log.(nodes), log.(min.(lowering, raising) ./ nodes))
        net_rates = rungs .* net.(log.(rungs))
        slower_rates = rungs .* exp.(slower.(log.(rungs)))
        lowering = slower_rates .+ max.(.-net_rates, 0)
        raising = slower_rates .+ max.(net_rates, 0)
    end
    PhononRateEquation(
        (lowering .+ rungs .* heating) .* time_unit^-1,
        (raising .+ rungs .* heating) .* time_unit^-1,
        time_unit,
    )
end

# The plain-number data of one mode from which the rungs are built, per jump operator:
# the recoil second moment α η₀² and the kicks η₀ h with their weights (none for
# dephasing, and none needed at first order).
struct RungLadder{T<:Unitful.Units}
    num_states::Int
    frame::Matrix{ComplexF64}
    beam_couplings::Vector{Matrix{ComplexF64}}
    lamb_dicke_factors::Vector{Float64}
    jump_operators::Vector{SparseMatrixCSC{ComplexF64,Int}}
    recoil_second_moments::Vector{Float64}
    recoil_kicks::Vector{Vector{Float64}}
    recoil_weights::Vector{Vector{Float64}}
    dissipator::SparseMatrixCSC{ComplexF64,Int}
    recoil_dissipator::SparseMatrixCSC{ComplexF64,Int}
    mode_frequency::Float64
    first_order::Bool
    bandwidth::Int
    time_unit::T
end

function RungLadder(model, mc; lamb_dicke_order, bandwidth, time_unit)
    first_order = lamb_dicke_order == 1
    !first_order &&
        mc.emission == :constant &&
        throw(
            ArgumentError(
                "Beyond first order in the Lamb–Dicke parameters the recoil kicks " *
                "need the emission-direction distribution; build the " *
                "MotionalCoupling with recoil_moment = :exact or :isotropic",
            ),
        )
    rate_matrix(M) = ComplexF64.(ustrip.(time_unit^-1, M))
    jump_operators = [
        sparse(ComplexF64.(ustrip.(time_unit^(-1 // 2), L))) for
        L in model.jump_operators
    ]
    recoil_second_moments = zeros(length(jump_operators))
    recoil_kicks = [Float64[] for _ in jump_operators]
    recoil_weights = [Float64[] for _ in jump_operators]
    for (i, k) in enumerate(mc.recoil_indices)
        η0 = mc.recoil_lamb_dicke[i]
        recoil_second_moments[k] = mc.recoil_moments[i] * η0^2
        if !first_order
            label = model.jump_labels[k]
            isotropic = mc.emission == :isotropic
            rank = isotropic ? 0 : multipole_rank(label.lower => label.upper)
            q = isotropic ? 0 : label.q
            projections, recoil_weights[k] = emission_rule(rank, q, mc.mode.direction)
            recoil_kicks[k] = η0 .* projections
        end
    end
    num_states = length(model.basis)
    RungLadder(
        num_states,
        Matrix{ComplexF64}(Diagonal(diag(rate_matrix(model.hamiltonian)))),
        rate_matrix.(model.couplings),
        mc.projected_lamb_dicke,
        jump_operators,
        recoil_second_moments,
        recoil_kicks,
        recoil_weights,
        sum(dissipator_superoperator, jump_operators),
        sum(dissipator_superoperator.(sqrt.(recoil_second_moments) .* jump_operators)),
        ustrip(time_unit^-1, mc.mode.frequency),
        first_order,
        bandwidth,
        time_unit,
    )
end

# Superoperators on column-stacked density matrices.
function commutator_superoperator(H)
    state_identity = sparse(ComplexF64(1) * I, size(H)...)
    kron(state_identity, sparse(H)) - kron(sparse(transpose(H)), state_identity)
end

function dissipator_superoperator(c)
    state_identity = sparse(ComplexF64(1) * I, size(c)...)
    decay = c' * c
    kron(conj(c), c) -
    (kron(state_identity, decay) + kron(sparse(transpose(decay)), state_identity)) / 2
end

# Nodes spaced by at most `node_ratio`, refined where a Fock coupling changes sign.
function rung_nodes(ladder::RungLadder, top_rung, node_ratio)
    nodes = [1]
    while last(nodes) < top_rung
        step = max(last(nodes) + 1, ceil(Int, node_ratio * last(nodes)))
        push!(nodes, min(top_rung, step))
    end
    if !ladder.first_order
        for η in ladder.lamb_dicke_factors, n in fock_coupling_zeros(η, top_rung)
            append!(nodes, clamp.(n .+ (-8, -4, -2, -1, 0, 1, 2, 4, 8), 1, top_rung))
        end
    end
    sort!(unique!(nodes))
end

# The Fock states after which the carrier or first-sideband element changes sign.
function fock_coupling_zeros(η, top_rung)
    changes_sign(k, n) =
        real(displacement_element(η, n, k) / im^k) *
        real(displacement_element(η, n + 1, k) / im^k) < 0
    [n for k in 0:1 for n in max(k, 1):(top_rung-1) if changes_sign(k, n)]
end

# The motional element ⟨n - k| e^{iη(a + a†)} |n⟩ at the order of the ladder.
function fock_coupling(ladder::RungLadder, η, n, k)
    ladder.first_order || return displacement_element(η, n, k)
    k == 0 ? 1.0 + 0im : k == 1 ? im * η * sqrt(n) : 0.0im
end

# The element repeated along the ladder of rung n: the carrier averaged over the rung.
function ladder_coupling(ladder::RungLadder, η, n, k)
    k == 0 &&
        return (fock_coupling(ladder, η, n, 0) + fock_coupling(ladder, η, n - 1, 0)) / 2
    k <= n ? fock_coupling(ladder, η, n, k) : 0.0im
end

# The terms H_k, k = 0 … bandwidth, of the ladder Hamiltonian H_0 ⊗ 1 + ∑_k H_k ⊗
# (E_-^k + E_+^k), with E_± the Fock shift operators.
function ladder_hamiltonians(ladder::RungLadder, n)
    map(0:ladder.bandwidth) do k
        H = k == 0 ? copy(ladder.frame) : zero(ladder.frame)
        for (C, η) in zip(ladder.beam_couplings, ladder.lamb_dicke_factors)
            element = ladder_coupling(ladder, η, n, k)
            H .+= (element .* C .+ conj(element) .* C') ./ 2
        end
        H
    end
end

# The probability per emission of changing the Fock state across rung n.
function recoil_fraction(ladder::RungLadder, n, jump)
    ladder.first_order && return n * ladder.recoil_second_moments[jump]
    kicks = abs2.(displacement_element.(ladder.recoil_kicks[jump], n, 1))
    dot(ladder.recoil_weights[jump], kicks)
end

# Returns the rates (Γ↓, Γ↑) of rung n in inverse time units: the nearest-neighbour
# steps with the net rate and the diffusion of the ladder, Γ↑/Γ↓ = (B + D)/(B - D).
# The window of four Fock states gives the ratio directly instead, to a few percent,
# which is the better source while the steps are directional (Γ↑/Γ↓ ≲ 0.8, i.e. B ≲
# 10|D|: the diffusion gives a small ratio as a small difference with ~1 % errors,
# which the steady occupation amplifies) and the window resolves the steps (a sideband
# reach 2|η|√n below one, the argument of the Bessel functions the displacement
# elements tend to). Nearer balance, the few percent in the window's ratio are the
# larger error, and the diffusion's ratio is exact at the balance point. Sub-Poissonian
# steps (a diffusion below the net rate) have no nearest-neighbour representation and
# take the window's ratio regardless. A window that disagrees with the ladder on the
# direction (close to where both change sign along the ladder) defers to the diffusion.
function rung_rates(ladder::RungLadder, n)
    (; net_rate, diffusion) = ladder_cumulants(ladder, n)
    total = diffusion
    directional = diffusion <= 10 * abs(net_rate)
    if diffusion < abs(net_rate) || (directional && sideband_reach(ladder, n) <= 1)
        ratio = stationary_rung_ratio(fock_window(ladder, n))
        if (ratio - 1) * net_rate > 0
            total = abs(net_rate) * (ratio + 1) / abs(ratio - 1)
        end
    end
    total = max(total, (1 + 1e-9) * abs(net_rate))
    ((total - net_rate) / 2, (total + net_rate) / 2)
end

# The sideband reach 2|η|√n of rung n, for the largest Lamb–Dicke factor of a beam or
# of a two-photon transition between two beams (whose sidebands reach as far).
function sideband_reach(ladder::RungLadder, n)
    η = ladder.lamb_dicke_factors
    largest = maximum(abs, η; init=0.0)
    for i in eachindex(η), j in (i+1):lastindex(η)
        largest = max(largest, abs(η[i] - η[j]))
    end
    2 * largest * sqrt(n)
end

ladder_net_rate(ladder::RungLadder, n) = ladder_cumulants(ladder, n).net_rate

# The first two cumulant rates of the phonon-number current of a ladder with the
# couplings of rung n throughout: `net_rate` Γ↑ - Γ↓ and `diffusion` Γ↑ + Γ↓. The
# sums over the Fock ladder ρ⁽ᵈ⁾ = ∑_m ⟨m + d|ρ|m⟩ close among themselves, for
# coherence orders d = -K … K. Counting the phonon number tilts the generator by
# e^{iχs} per term, with s the change of m + d/2 it causes, L(χ) = L + iχ C₁ -
# χ² C₂/2 + …, so that with the stationary ρ and the traceless x solving L x = -(C₁
# - D) ρ, the cumulants are D = tr C₁ρ and B = tr C₂ρ + 2 tr C₁x. The first-order
# ladder carries the recoil as the jumps c ⊗ η₀√n(E₊ + E₋), counted likewise; the
# full-order ladder does not, and the recoil adds its diffusion 2n ∑ m₂ R per jump.
function ladder_cumulants(ladder::RungLadder, n)
    hamiltonians = ladder_hamiltonians(ladder, n)
    num_states = ladder.num_states
    K = ladder.bandwidth
    num_orders = 2K + 1
    vector_size = num_states^2
    same_order = sparse(1.0I, num_orders, num_orders)
    from_above(k) = spdiagm(k => ones(num_orders - k))   # ρ⁽ᵈ⁾ from ρ⁽ᵈ⁺ᵏ⁾
    from_below(k) = spdiagm(-k => ones(num_orders - k))  # ρ⁽ᵈ⁾ from ρ⁽ᵈ⁻ᵏ⁾
    state_identity = sparse(ComplexF64(1) * I, num_states, num_states)
    L =
        kron(
            same_order,
            ladder.dissipator - im * commutator_superoperator(hamiltonians[1]),
        ) + kron(
            spdiagm(0 => -im * ladder.mode_frequency .* ((-K):K)),
            sparse(I, vector_size, vector_size),
        )
    C1 = spzeros(ComplexF64, size(L)...)
    C2 = spzeros(ComplexF64, size(L)...)
    for k in 1:K
        H = sparse(hamiltonians[k+1])
        left = -im * kron(state_identity, H)                     # -i H_k ρ
        right = im * kron(sparse(transpose(H)), state_identity)  # +i ρ H_k
        L += kron(from_above(k) + from_below(k), left + right)
        # E₋ᵏ on the ket or E₊ᵏ on the bra come from order d + k and lower the
        # counted number by k/2; their mirror images raise it.
        C1 +=
            (k / 2) *
            (kron(from_above(k), right - left) + kron(from_below(k), left - right))
        C2 += (k^2 / 4) * kron(from_above(k) + from_below(k), left + right)
    end
    if ladder.first_order
        shift = from_above(2) + from_below(2)
        L += n * kron(2 * same_order + shift, ladder.recoil_dissipator)
        # The kicks c ρ c† move both indices by ±1 (counted ±1, same order) or
        # opposite ways (uncounted, order ± 2); the decay terms with E₊² or E₋² on
        # one side count ±1 and change the order by 2.
        kicks = sum(
            m2 * kron(conj(c), c) for
            (c, m2) in zip(ladder.jump_operators, ladder.recoil_second_moments)
        )
        decay = sum(
            m2 * c' * c for
            (c, m2) in zip(ladder.jump_operators, ladder.recoil_second_moments)
        )
        decay_left = kron(state_identity, decay)
        decay_right = kron(sparse(transpose(decay)), state_identity)
        C1 -= (n / 2) * kron(from_below(2) - from_above(2), decay_left - decay_right)
        C2 +=
            n *
            (2 * kron(same_order, kicks) - kron(shift, decay_left + decay_right) / 2)
    end
    order_range(d) = (1:vector_size) .+ (d + K) * vector_size
    diagonal = collect(order_range(0)[1:(num_states+1):end])
    (; F, row) = trace_row_system(L, diagonal)
    unit_trace = zeros(ComplexF64, size(L, 1))
    unit_trace[row] = 1
    ρ = F \ unit_trace
    trace(v) = sum(v[diagonal])
    net_rate = real(trace(C1 * ρ))
    traceless = -(C1 * ρ .- net_rate .* ρ)
    traceless[row] = 0
    x = F \ traceless
    finalize(F)
    diffusion = real(trace(C2 * ρ) + 2 * trace(C1 * x))
    if !ladder.first_order
        ρ0 = reshape(ρ[order_range(0)], num_states, num_states)
        diffusion += sum(
            real(tr(c' * c * ρ0)) * 2n * m2 for
            (c, m2) in zip(ladder.jump_operators, ladder.recoil_second_moments)
        )
    end
    (; net_rate, diffusion)
end

# The null vector of the generator L whose entries `diagonal` sum to one.
function stationary_vector(L, diagonal)
    (; F, row) = trace_row_system(L, diagonal)
    unit_trace = zeros(ComplexF64, size(L, 1))
    unit_trace[row] = 1
    vector = F \ unit_trace
    finalize(F)
    vector
end

# The closed internal ⊗ Fock master equation of the states n - 2 … n + 1 around rung n.
# A kick between Fock states keeps the total emission rate: its probability is taken
# from the unkicked jump (except in the first-order model, whose rate grows with n).
function fock_window(ladder::RungLadder, n)
    num_states = ladder.num_states
    fock_states = max(0, n-2):(n+1)
    num_fock = length(fock_states)
    fock_identity = sparse(ComplexF64(1) * I, num_fock, num_fock)
    fock_energies = spdiagm(0 => ComplexF64.(fock_states .- first(fock_states)))
    H =
        kron(sparse(ladder.frame), fock_identity) +
        ladder.mode_frequency *
        kron(sparse(ComplexF64(1) * I, num_states, num_states), fock_energies)
    for (C, η) in zip(ladder.beam_couplings, ladder.lamb_dicke_factors)
        motional = [
            fock_coupling(ladder, η, max(l, m), abs(l - m)) for
            l in fock_states, m in fock_states
        ]
        coupling = kron(sparse(C), sparse(motional)) / 2
        H += coupling + coupling'
    end
    jumps = SparseMatrixCSC{ComplexF64,Int}[]
    for (jump, c) in enumerate(ladder.jump_operators)
        kicks = [recoil_fraction(ladder, m, jump) for m in fock_states[2:end]]
        unkicked =
            ladder.first_order ? ones(num_fock) :
            max.(0.0, 1 .- [kicks; 0.0] .- [0.0; kicks])
        push!(jumps, kron(c, spdiagm(0 => sqrt.(unkicked))))
        any(>(0), kicks) &&
            push!(jumps, kron(c, spdiagm(1 => sqrt.(kicks), -1 => sqrt.(kicks))))
    end
    (;
        liouvillian=BandLiouvillian(
            H,
            jumps,
            num_states,
            num_fock;
            bandwidth=ladder.bandwidth,
            time_unit=ladder.time_unit,
        ),
        upper_position=n - first(fock_states) + 1,
    )
end

function stationary_rung_ratio(window)
    bl = window.liouvillian
    steady = stationary_vector(bl.L, [bl.index[r, r] for r in 1:full_dimension(bl)])
    populations = fock_populations(bl, steady, 1)
    populations[window.upper_position] / populations[window.upper_position-1]
end

phonon_numbers(r::PhononRateEquation) = 0:length(r.lowering_rates)

initial_populations(r::PhononRateEquation, initial_nbar::Real) =
    thermal_populations(initial_nbar, length(r.lowering_rates) + 1)
initial_populations(::PhononRateEquation, populations::AbstractVector) =
    collect(Float64, populations)

# The generator of dp/dt in inverse time units.
function rate_generator(r::PhononRateEquation)
    lowering = ustrip.(r.time_unit^-1, r.lowering_rates)
    raising = ustrip.(r.time_unit^-1, r.raising_rates)
    Tridiagonal(raising, -[raising; 0] .- [0; lowering], lowering)
end

"""
    fock_populations(r::PhononRateEquation) -> Vector{Float64}

Returns the stationary phonon populations, by detailed balance across each rung.
"""
function fock_populations(r::PhononRateEquation)
    populations = exp.(cumsum([0.0; log.(r.raising_rates ./ r.lowering_rates)]))
    populations ./ sum(populations)
end

"""
    steady_value(r::PhononRateEquation) -> Float64

Returns the stationary mean phonon number.
"""
steady_value(r::PhononRateEquation) = dot(phonon_numbers(r), fock_populations(r))

"""
    integral_relaxation_time(r::PhononRateEquation, initial) -> Unitful.Time

Returns the integral relaxation time of the mean phonon number, as for
[`IntegratedTransient`](@ref), from `initial`: a thermal mean occupation or a
vector of phonon populations.
"""
function integral_relaxation_time(r::PhononRateEquation, initial)
    start = initial_populations(r, initial)
    steady = fock_populations(r)
    system = Matrix(rate_generator(r))
    system[1, :] .= 1   # the transient integral is traceless
    transient = steady .- start
    transient[1] = 0
    area = dot(phonon_numbers(r), system \ transient)
    area / (dot(phonon_numbers(r), start) - steady_value(r)) * r.time_unit
end

function rate_problem(r::PhononRateEquation, initial, time_span)
    generator = rate_generator(r)
    ODEProblem(
        (dp, p, _, _) -> mul!(dp, generator, p),
        initial_populations(r, initial),
        time_span,
    )
end

"""
    cooling_curve(r::PhononRateEquation, initial, ts; fock = false, reltol = 1e-8, abstol = 1e-10, alg = Vern7(), kwargs...)

Integrates the rate equation from `initial` (a thermal mean occupation or a
vector of phonon populations) and returns `(; t, nbar, fock)` as for a
[`BandLiouvillian`](@ref), with one column in `nbar` and one matrix in `fock`.
"""
function cooling_curve(
    r::PhononRateEquation,
    initial,
    ts;
    fock::Bool=false,
    reltol=1e-8,
    abstol=1e-10,
    alg=Vern7(),
    kwargs...,
)
    times = eltype(ts) <: Quantity ? ustrip.(r.time_unit, ts) : collect(float.(ts))
    problem = rate_problem(r, initial, (first(times), last(times)))
    solution = solve(problem, alg; saveat=times, reltol, abstol, kwargs...)
    successful_retcode(solution) ||
        error("Rate-equation propagation failed: $(solution.retcode)")
    populations = reduce(hcat, solution.u)
    (;
        t=collect(solution.t) .* r.time_unit,
        nbar=reshape(phonon_numbers(r)' * populations, :, 1),
        fock=fock ? [populations] : nothing,
    )
end

"""
    time_to_nbar(r::PhononRateEquation, initial, target_nbar; reltol = 1e-8, abstol = 1e-10, alg = Vern7()) -> Unitful.Time

Returns the time at which the mean phonon number first reaches `target_nbar`
when cooling from `initial` (a thermal mean occupation or a vector of phonon
populations): zero if it starts there or below, `Inf` if its steady value lies
above.
"""
function time_to_nbar(
    r::PhononRateEquation,
    initial,
    target_nbar::Real;
    reltol=1e-8,
    abstol=1e-10,
    alg=Vern7(),
)
    start = initial_populations(r, initial)
    dot(phonon_numbers(r), start) <= target_nbar && return 0.0 * r.time_unit
    steady_value(r) >= target_nbar && return Inf * r.time_unit
    reached = ContinuousCallback(
        (p, _, _) -> dot(phonon_numbers(r), p) - target_nbar,
        terminate!,
    )
    solution = solve(
        rate_problem(r, start, (0.0, Inf)),
        alg;
        callback=reached,
        save_everystep=false,
        reltol,
        abstol,
    )
    successful_retcode(solution) ||
        error("Rate-equation propagation failed: $(solution.retcode)")
    last(solution.t) * r.time_unit
end

function cooling_metrics(
    model::LindbladModel,
    mcs::AbstractVector{<:MotionalCoupling},
    method::SaturatedRates;
    heating_rate=0.0u"s^-1",
    ρ=nothing,
    time_unit::Unitful.Units=u"µs",
)
    num_modes = length(mcs)
    per_mode(x) = x isa Real ? fill(x, num_modes) : x
    initial_nbars = per_mode(method.initial_nbar)
    num_focks =
        isnothing(method.num_fock) ? fock_truncation.(initial_nbars) :
        per_mode(method.num_fock)
    length(initial_nbars) == length(num_focks) == num_modes || throw(
        ArgumentError("initial_nbar and num_fock must be one value or one per mode"),
    )
    map(
        mcs,
        initial_nbars,
        num_focks,
        mode_heating_rates(heating_rate, num_modes),
    ) do mc, initial_nbar, num_fock, rate
        equation = PhononRateEquation(
            model,
            mc;
            num_fock,
            lamb_dicke_order=method.lamb_dicke_order,
            bandwidth=method.bandwidth,
            node_ratio=method.node_ratio,
            heating_rate=rate,
            time_unit,
        )
        CoolingMetrics(
            method,
            steady_value(equation),
            integral_relaxation_time(equation, initial_nbar),
            nothing,
            nothing,
            nothing,
            equation,
        )
    end
end

export SaturatedRates, PhononRateEquation, time_to_nbar
