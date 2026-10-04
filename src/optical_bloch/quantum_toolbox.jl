# QuantumToolbox.jl integration: the function stubs (and user-facing
# docstrings) of the solvers implemented in the LevelsQuantumToolboxExt
# package extension on top of the LindbladModel/MotionalCoupling data. The
# extension also adds methods to QuantumToolbox's own `liouvillian` and
# `steadystate` for a LindbladModel (documented there).

"""
    populations(ρ, model::LindbladModel) -> Vector{Float64}
    populations(ρ, model::LindbladModel, level) -> Float64

Returns the populations of the basis states in the (internal or reduced
internal) density matrix `ρ`, or the total population of the given level (a
fine-structure level covering all its hyperfine levels, cf.
[`level_projector`](@ref)).

!!! note
    Defined in the `LevelsQuantumToolboxExt` package extension, i.e. only once
    QuantumToolbox.jl is loaded.
"""
function populations end

"""
    cooling_rates(model::LindbladModel, mc::MotionalCoupling; ρ = steadystate(model), heating_rate = 0)
    cooling_rates(model::LindbladModel, mcs::AbstractVector{<:MotionalCoupling}; ρ, heating_rate)

Returns the laser-cooling rate coefficients of the mode(s) by adiabatic
elimination of the internal dynamics — valid when the internal relaxation is
fast compared to the motional dynamics (``η Ω ≪ Γ, ω_m``), as a `NamedTuple`
`(; A_plus, A_minus, nbar, τ_c)` (a vector of them for several modes):

```math
A_± = 2 \\, \\mathrm{Re} \\, S(∓ω_m) + D, \\qquad
S(ω) = \\mathrm{Tr}\\left[ H_\\mathrm{sb} \\, (-\\mathcal{L} + i ω)^{-1} \\left( H_\\mathrm{sb} ρ_∞ \\right) \\right],
\\qquad
D = \\sum_k α_k η_{0,k}^2 \\, \\mathrm{Tr}\\left[ L_k^† L_k ρ_∞ \\right],
```

with ``\\mathcal{L}`` the internal Liouvillian, ``ρ_∞`` its steady state,
``H_\\mathrm{sb}`` the sideband Hamiltonian and ``D`` the recoil diffusion from
spontaneous emission (cf. [`MotionalCoupling`](@ref); Cirac et al., Phys. Rev.
A **46**, 2668 (1992); Morigi, Phys. Rev. A **67**, 033402 (2003)). ``A_+``
raises and ``A_-`` lowers the phonon number, so the mean occupation relaxes as
``\\dot{\\bar n} = A_+ (\\bar n + 1) - A_- \\bar n`` (plus an optional
`heating_rate` in quanta per time) towards
``\\bar n_∞ = (A_+ + \\dot{n}_\\mathrm{heat}) / (A_- - A_+)`` with the time
constant ``τ_c = 1 / (A_- - A_+)``; both are `NaN` when the mode is heated
(``A_- ≤ A_+``). The rates are in inverse `time_unit` (µs⁻¹ by default), `nbar`
dimensionless, `τ_c` in `time_unit`.

For several modes the internal steady state and Liouvillian factorisation are
shared; the modes are independent at this order.

!!! note
    Defined in the `LevelsQuantumToolboxExt` package extension, i.e. only once
    QuantumToolbox.jl is loaded.
"""
function cooling_rates end

"""
    motional_model(model::LindbladModel, mc::MotionalCoupling; num_fock, time_unit = u"µs")
    motional_model(model::LindbladModel, mcs::AbstractVector{<:MotionalCoupling}; num_fock, time_unit)

Returns the full internal ⊗ motional model as QuantumToolbox operators,
`(; H, c_ops, dims)`: the Hamiltonian
``H_\\mathrm{int} ⊗ 1 + \\sum_m ω_m a_m^† a_m + \\sum_m H_{\\mathrm{sb},m} ⊗ (a_m + a_m^†)``
on the internal states and one Fock space per mode truncated at `num_fock` levels
(a number, or one per mode), and the jump operators ``L ⊗ 1`` together with the
recoil operators ``\\sqrt{α} η_0 L ⊗ (a_m + a_m^†)`` of every mode. Recoil kicks
are treated independently per mode (the emission-pattern cross moments between
two modes are dropped). Requires a static internal model (no beat harmonics).

!!! note
    Defined in the `LevelsQuantumToolboxExt` package extension, i.e. only once
    QuantumToolbox.jl is loaded.
"""
function motional_model end

"""
    mean_phonon_number(ρ, model::LindbladModel, mcs; method = :direct)

Returns the mean phonon number of each mode from a density matrix of the full
[`motional_model`](@ref) (a vector for several modes, a number for one), from
the populations ``p_n`` of the reduced motional state, with `method`

  - `:direct`: the expectation value ``⟨a^† a⟩ = \\sum_n n p_n``, the actual
    observable, but biased low unless the truncation is loose
    (``p_{n_\\mathrm{max}-1} ≪ 1``);
  - `:thermal_ratio`: ``\\bar n = r/(1 - r)`` from the ratio estimator
    ``r = \\sum_{n ≥ 1} p_n / \\sum_{n ≤ n_\\mathrm{max}-2} p_n`` of a thermal
    distribution ``p_n ∝ r^n``, which is exact at any truncation if the state
    is thermal. That is the case in the rate-equation limit, where the phonon
    populations form a birth–death chain in detailed balance, whose steady
    state a truncation merely renormalises. Beyond it (``η Ω`` not ≪ ``ω_m``),
    ``p_{n+1}/p_n`` grows with ``n``, and for a loose truncation the estimate
    tends to ``1/p_0 - 1`` rather than ``⟨a^† a⟩`` (≈ 10 % low for ⁴⁰Ca⁺ Λ
    dark-resonance cooling, from ``η Ω ≈ ω_m/4`` up). It is meant for
    truncations too tight for `:direct` to converge.

!!! note
    Defined in the `LevelsQuantumToolboxExt` package extension, i.e. only once
    QuantumToolbox.jl is loaded.
"""
function mean_phonon_number end

"""
    cooling_time(H, c_ops; time_unit = u"µs", eigvals = 12, dense_limit = 2500, weight_threshold = 0.1)

Returns the characteristic time of the slowest relaxation of the mean phonon
number under the master equation with Hamiltonian `H` and jump operators
`c_ops` (QuantumToolbox operators on an internal ⊗ motional space, e.g. from
[`motional_model`](@ref)): the negative inverse real part of the slowest
Liouvillian eigenvalue, after the stationary state, whose eigenmode has
a weight ``|\\mathrm{Tr}[\\hat N v]| / ‖v‖`` in the total phonon number of at
least `weight_threshold` times the largest such weight. Coherence modes — which
in the rate-equation limit relax at *half* the rate of ``\\bar n`` and would
otherwise be picked up as the slowest, carrying only a weight of order ``η²``
through the sideband coupling — are thereby skipped, so the result matches the
`τ_c` of [`cooling_rates`](@ref) where the adiabatic elimination holds (and the
Fock truncation is loose enough not to shape the slowest mode). Liouvillians of dimension up to `dense_limit` are
diagonalised densely, larger ones with a shift-inverted Arnoldi iteration for
`eigvals` eigenpairs around zero.

!!! note
    Defined in the `LevelsQuantumToolboxExt` package extension, i.e. only once
    QuantumToolbox.jl is loaded.
"""
function cooling_time end

export populations, cooling_rates, motional_model, mean_phonon_number, cooling_time
