# The motional layer: Lamb–Dicke factors, the first-order sideband Hamiltonian
# and the recoil jump operators of one motional mode, derived from a
# LindbladModel of the internal states.

"""
Coupling of one [`MotionalMode`](@ref) to the internal dynamics of a
[`LindbladModel`](@ref), to first order in the Lamb–Dicke parameters.

With ``x = x_0 (a + a^†)``, ``x_0 = \\sqrt{ħ / 2 M ω_m}``, along the mode
direction, the full Hamiltonian is ``H_\\mathrm{int} ⊗ 1 + ω_m a^† a +
H_\\mathrm{sb} ⊗ (a + a^†)`` and each spontaneous-emission jump operator ``L``
acquires the companion ``\\sqrt{α} η_0 L ⊗ (a + a^†)`` for the recoil of the
emitted photon.

Two distinct Lamb–Dicke quantities appear, and the field names spell out which:

- `projected_lamb_dicke[l]` ``= (\\vec{k}_l ⋅ \\vec{e}_m) x_0`` is the **signed,
  projected** Lamb–Dicke factor of beam ``l`` — the coefficient of ``(a + a^†)``
  in ``e^{i \\vec{k}_l ⋅ \\vec{r}}`` — which is what enters the sideband
  Hamiltonian ``H_\\mathrm{sb} = \\sum_l η_l \\, i (C_l - C_l^†)/2``. For a Raman
  process through two beams, the difference of their values is the two-photon
  factor.
- `recoil_lamb_dicke[k]` ``= |\\vec{k}| x_0`` is the **unprojected** Lamb–Dicke
  parameter of the decay transition behind jump operator ``k``; the direction of
  the emitted photon is averaged over its radiation pattern into the second
  moment ``α_k = ⟨(\\hat{k} ⋅ \\vec{e}_m)^2⟩`` (`recoil_moments`), which depends
  on the emitted spherical component ``q`` and the angle between the mode and
  the quantisation axis: ``1/5`` (``q = 0``) and ``2/5`` (``|q| = 1``) for a mode
  along ẑ, ``1/3`` on average over ``q`` (cf. [`recoil_moment`](@ref)).

`recoil_operators[k]` is ``\\sqrt{α_k} η_{0,k} L_k`` for the ``k``-th
spontaneous-emission operator of the model (`recoil_indices[k]` its index in
`model.jump_operators`; laser-dephasing operators carry no recoil).
`sideband_harmonics` mirrors `model.harmonics` for beams with a beat note.

Constructed via [`motional_coupling`](@ref).
"""
struct MotionalCoupling{M<:MotionalMode,H<:Quantity,W<:Quantity,J<:Quantity}
    mode::M
    projected_lamb_dicke::Vector{Float64}
    recoil_lamb_dicke::Vector{Float64}
    recoil_moments::Vector{Float64}
    sideband_hamiltonian::Matrix{H}
    sideband_harmonics::Vector{Tuple{W,Matrix{H}}}
    recoil_operators::Vector{Matrix{J}}
    recoil_indices::Vector{Int}
end

"""
    lamb_dicke(species, mode::MotionalMode, wavevector; projected = true)
    lamb_dicke(species, mode::MotionalMode, beam::LaserBeam; projected = true)

Returns the Lamb–Dicke factor of a plane wave with the given wavevector (a
Cartesian 3-vector of inverse lengths, or the beam's ``k = ω n / c``) and the
given mode of an ion of the species' mass: the signed **projected** value
``(\\vec{k} ⋅ \\vec{e}_m) \\sqrt{ħ / 2 M ω_m}`` by default, or with
`projected = false` the **unprojected** ``|\\vec{k}| \\sqrt{ħ / 2 M ω_m}``, the
conventional Lamb–Dicke parameter of the transition (cf.
[`MotionalCoupling`](@ref) for which is which).
"""
function lamb_dicke(species, mode::MotionalMode, wavevector; projected::Bool=true)
    x0 = sqrt(u"ħ" / (2 * species.mass * mode.frequency))
    k = projected ? sum(wavevector .* mode.direction) : sqrt(sum(abs2, wavevector))
    uconvert(NoUnits, k * x0)
end

function lamb_dicke(species, mode::MotionalMode, beam::LaserBeam; projected::Bool=true)
    k = photon_energy(species, beam.frequency) / (u"ħ" * u"c")
    lamb_dicke(species, mode, k .* beam.n; projected)
end

# Angular radiation patterns of the electric multipole components, as
# polynomial coefficients in c = cos θ (θ from the quantisation axis), up to a
# constant: E1 q = 0: sin²θ, |q| = 1: (1 + cos²θ)/2; E2 q = 0: sin²θ cos²θ,
# |q| = 1: 1 − 3 cos²θ + 4 cos⁴θ, |q| = 2: 1 − cos⁴θ (Jackson, Classical
# Electrodynamics, Sec. 9.7).
function radiation_pattern(rank, q)
    if rank == 1
        q == 0 ? (1.0, 0.0, -1.0) : (1.0, 0.0, 1.0)
    elseif rank == 2
        if q == 0
            (0.0, 0.0, 1.0, 0.0, -1.0)
        elseif abs(q) == 1
            (1.0, 0.0, -3.0, 0.0, 4.0)
        else
            (1.0, 0.0, 0.0, 0.0, -1.0)
        end
    else
        throw(ArgumentError("Radiation patterns are tabulated for ranks 1 and 2 only"))
    end
end

# ∫_{-1}^{1} c^n p(c) dc for a polynomial p given by its coefficients.
function moment_integral(coefficients, n)
    total = 0.0
    for (i, a) in enumerate(coefficients)
        power = n + i - 1
        iseven(power) && (total += a * 2 / (power + 1))
    end
    total
end

"""
    recoil_moment(rank, q, direction)

Returns the second moment ``α = ⟨(\\hat{k} ⋅ \\vec{e}_m)^2⟩`` of the direction of a
photon emitted in the spherical component `q` of an electric multipole
transition of the given `rank`, projected onto the unit mode `direction`
(quantisation axis along ẑ): the fraction of the squared recoil that a
spontaneous emission imparts, on average, to that mode.

By the azimuthal symmetry of the radiation patterns,
``α = ⟨\\hat{k}_z^2⟩ \\cos^2 θ_m + \\tfrac{1}{2}(1 - ⟨\\hat{k}_z^2⟩) \\sin^2 θ_m``
with ``θ_m`` the angle between mode and quantisation axis; for electric-dipole
emission ``⟨\\hat{k}_z^2⟩ = 1/5`` (``q = 0``) and ``2/5`` (``|q| = 1``), so a mode
along ẑ sees ``1/5`` and ``2/5``, one transverse to it ``2/5`` and ``3/10``,
and the average over the three components is the isotropic ``1/3`` in every
direction (cf. Javanainen & Stenholm, Appl. Phys. **21**, 35 (1980)).
"""
function recoil_moment(rank, q, direction)
    abs(q) <= rank || throw(ArgumentError("|q| must not exceed the rank"))
    p = radiation_pattern(rank, q)
    kz2 = moment_integral(p, 2) / moment_integral(p, 0)
    cos2 = abs2(direction[3]) / sum(abs2, direction)
    kz2 * cos2 + (1 - kz2) / 2 * (1 - cos2)
end

# The recoil second moment of one decay operator per the `recoil_moment`
# setting: exact from the radiation pattern of its emitted component, or a
# constant.
decay_recoil_moment(setting::Symbol, label::DecayLabel, mode::MotionalMode) =
    setting === :exact ?
    recoil_moment(multipole_rank(label.lower, label.upper), label.q, mode.direction) :
    throw(ArgumentError("recoil_moment must be :exact or a number, got $setting"))
decay_recoil_moment(setting::Real, ::DecayLabel, ::MotionalMode) = Float64(setting)

"""
    motional_coupling(species, scheme::LaserScheme, model::LindbladModel, mode; recoil_moment = :exact)
    motional_coupling(species, scheme, model, modes::AbstractVector) -> Vector{MotionalCoupling}

Derives the [`MotionalCoupling`](@ref) of one [`MotionalMode`](@ref) (or of each
of several) to the internal dynamics of `model`, built from `scheme`: the
projected Lamb–Dicke factor of each beam, the sideband Hamiltonian
``\\sum_l η_l \\, i (C_l - C_l^†)/2`` (its beat-note harmonics kept alongside),
and the recoil operators ``\\sqrt{α} η_0 L`` of the spontaneous-emission jump
operators. `recoil_moment = :exact` evaluates ``α`` per emitted component from
the radiation patterns ([`recoil_moment`](@ref)); a number uses that constant
for every decay (e.g. the often-quoted ``2/5``).
"""
function motional_coupling(
    species,
    scheme::LaserScheme,
    model::LindbladModel,
    mode::MotionalMode;
    recoil_moment=:exact,
)
    basis = scheme.basis
    B = scheme.static_field
    n = length(basis)

    # Sideband Hamiltonian: coefficient of (a + a†) from e^{ik·r} ≈ 1 + iη(a + a†).
    η = [lamb_dicke(species, mode, beam) for beam in scheme.beams]
    H_sb = zeros(eltype(model.hamiltonian), n, n)
    beat_of = Dict(b => w for (b, w) in model.frame.beats)
    harmonics = Tuple{typeof(1.0 * ANGULAR_UNIT),typeof(H_sb)}[]
    for (b, beam) in enumerate(scheme.beams)
        lo, hi = beam_levels(beam)
        C = coupling_matrix(species, basis, lo => hi, beam.intensity, beam.ε, beam.n, B)
        term = (η[b] * im) .* (C .- C') ./ 2
        w = get(beat_of, b, zero(1.0 * ANGULAR_UNIT))
        if iszero(w)
            H_sb .+= term
        elseif w > zero(w)
            push!(harmonics, (w, (η[b] * im) .* C ./ 2))
        else
            push!(harmonics, (-w, (η[b] * im) .* Matrix(C') ./ 2))
        end
    end

    # Recoil: the unprojected Lamb–Dicke parameter of each decay transition,
    # weighted by the emission-pattern second moment along the mode.
    recoil_indices = findall(l -> l isa DecayLabel, model.jump_labels)
    η0 = Float64[]
    α = Float64[]
    ops = similar(model.jump_operators, 0)
    for k in recoil_indices
        label = model.jump_labels[k]
        k_photon = transition_frequency(species, label.lower, label.upper) / u"c"
        η0_k = lamb_dicke(species, mode, k_photon .* mode.direction; projected=false)
        α_k = decay_recoil_moment(recoil_moment, label, mode)
        push!(η0, η0_k)
        push!(α, α_k)
        push!(ops, (sqrt(α_k) * η0_k) .* model.jump_operators[k])
    end

    MotionalCoupling(mode, η, η0, α, H_sb, harmonics, ops, recoil_indices)
end

motional_coupling(
    species,
    scheme::LaserScheme,
    model::LindbladModel,
    modes::AbstractVector;
    recoil_moment=:exact,
) = [motional_coupling(species, scheme, model, mode; recoil_moment) for mode in modes]

export MotionalCoupling, lamb_dicke, recoil_moment, motional_coupling
