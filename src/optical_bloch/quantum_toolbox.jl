# QuantumToolbox.jl integration: the types and function stubs (with the
# user-facing docstrings) of the solvers implemented in the
# LevelsQuantumToolboxExt package extension on top of the
# LindbladModel/MotionalCoupling data. The extension also adds methods to
# QuantumToolbox's own `liouvillian` and `steadystate` for a LindbladModel and
# a MotionalModel (documented there).

"""
Label of a recoil jump operator of a [`MotionalModel`](@ref): the
[`DecayLabel`](@ref) of the spontaneous-emission operator it derives from, the
`mode` it kicks (first-order model: one operator ``\\sqrt{α} η_0 L ⊗ (a_m +
a_m^†)`` per mode, `node = 0`) or `mode = 0` for the joint displacement
operators ``\\sqrt{w_j} L ⊗ D(-η_0 \\vec{h}_j)`` of the emission-direction
quadrature (`node` ``= j``, cf. [`emission_rule`](@ref)) of the higher-order
models.
"""
struct RecoilLabel
    decay::DecayLabel
    mode::Int
    node::Int
end

"""
Label of a heating-bath jump operator ``\\sqrt{\\dot n} a_m`` (`raising =
false`) or ``\\sqrt{\\dot n} a_m^†`` (`raising = true`) of mode `mode` in a
[`MotionalModel`](@ref).
"""
struct HeatingLabel
    mode::Int
    raising::Bool
end

"""
The full internal ⊗ motional master equation of a [`LindbladModel`](@ref) and
the [`MotionalCoupling`](@ref)s of its modes, as built by
[`motional_model`](@ref): the Hamiltonian `H` and jump operators `c_ops`
(QuantumToolbox operators on the internal states ⊗ one Fock space per mode,
internal factor first, in inverse `time_unit`), one label per jump operator in
`c_labels` (the [`DecayLabel`](@ref) or [`DephasingLabel`](@ref) of an
internal operator acting as ``L ⊗ 1``, a [`RecoilLabel`](@ref) for the recoil
operators of a decay, a [`HeatingLabel`](@ref) for the bath operators), the
layout — `num_states` internal states, `num_fock` Fock states per mode, the
`modes` — and the `lamb_dicke_order` the model was expanded to.

The extension defines `steadystate(mm)`, `liouvillian(mm)`,
[`mean_phonon_number`](@ref), [`fock_populations`](@ref),
[`thermal_state`](@ref), [`cooling_time`](@ref), [`cooling_curve`](@ref),
[`BandLiouvillian`](@ref) and [`IntegratedTransientSolver`](@ref) on it.
"""
struct MotionalModel{O,C,M<:MotionalMode,T<:Unitful.Units}
    H::O
    c_ops::Vector{C}
    c_labels::Vector{Any}
    num_states::Int
    num_fock::Vector{Int}
    modes::Vector{M}
    lamb_dicke_order::Float64
    time_unit::T
end

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
    motional_model(model::LindbladModel, mc::MotionalCoupling; num_fock, lamb_dicke_order = 1, heating_rate = 0, nodes, time_unit = u"µs")
    motional_model(model::LindbladModel, mcs::AbstractVector{<:MotionalCoupling}; num_fock, …) -> MotionalModel

Returns the full internal ⊗ motional master equation as a
[`MotionalModel`](@ref): the internal states and one Fock space per mode
truncated at `num_fock` states (a number, or one per mode), with the
Hamiltonian

```math
H = H_\\mathrm{int} ⊗ 1 + \\sum_m ω_m a_m^† a_m
  + \\sum_b \\tfrac{1}{2} \\left( C_b ⊗ D_b + \\mathrm{h.c.} \\right), \\qquad
D_b = e^{i \\sum_m η_{b,m} (a_m + a_m^†)},
```

where ``C_b`` is the coupling matrix of beam ``b`` and ``η_{b,m}`` its
projected Lamb–Dicke factor on mode ``m`` (cf. [`MotionalCoupling`](@ref)) —
the plane-wave phase ``e^{i\\vec{k}_b ⋅ \\vec{r}}`` of the beam acting on the
motion — and, for every spontaneous-emission operator ``L``, the recoil of the
emitted photon, ``L ⊗ e^{-i \\vec{k} ⋅ \\vec{r}}`` averaged over the emission
direction distribution of its multipole component (Kulosa et al., New J. Phys.
**25**, 053008 (2023), Eqs. (6), (7), (A.3)).

`lamb_dicke_order` sets the Taylor order ``k`` to which these displacement
operators are expanded, `Inf` keeping them exactly
([`displacement_elements`](@ref)):

  - `1` (default) is the Lamb–Dicke-regime model of the laser-cooling
    literature, ``D_b ≈ 1 + i\\sum_m η_{b,m}(a_m + a_m^†)`` — the sideband
    Hamiltonian ``H_{\\mathrm{sb},m} ⊗ (a_m + a_m^†)`` of the coupling — and,
    per decay operator and mode, the recoil operator ``\\sqrt{α_m} η_{0,m} L ⊗
    (a_m + a_m^†)`` next to ``L ⊗ 1`` (Cirac et al., Phys. Rev. A **46**, 2668
    (1992); Morigi, Phys. Rev. A **67**, 033402 (2003)). It is the only order
    closed in the second moment ``α`` of the emission distribution alone
    (hence the only one accepting a constant `recoil_moment`), it reproduces
    the recoil heating ``α η_0^2`` per emission exactly, and it is the model
    whose adiabatic elimination gives [`cooling_rates`](@ref). Recoil
    correlations between modes are dropped at this order.
  - A finite ``k ≥ 2`` truncates ``D_b`` and the recoil displacements at total
    order ``k`` in the Lamb–Dicke factors (the elements of the truncated
    operators are exact), with the recoil averaged over an emission-direction
    quadrature ([`emission_rule`](@ref)) of ``k + 1`` nodes per decay (exact for
    the moments entering), which requires the [`MotionalCoupling`](@ref)s to
    carry the emission distribution (`recoil_moment = :exact` or
    `:isotropic`). Successive orders converge to the exact model in powers of
    ``η\\sqrt{n}`` — the recoil heating per emission, exact at first order, is
    recovered only to ``O(η^{k+2})`` — which is why none of them is cheaper
    where it is accurate.
  - `Inf` keeps the full operators: beam couplings through the displacement
    matrix elements, each decay operator replaced by the quadrature set
    ``\\sqrt{w_j} L ⊗ D(-η_0 \\vec{h}_j)`` over `nodes` emission directions (12
    by default; a product rule on the sphere for several modes, which makes
    the cross-mode recoil correlations exact). Required for
    ``η\\sqrt{2\\bar n + 1} ≳ 0.5`` — a Doppler-cooled ion at ``η ≈ 0.1`` — where
    the first-order model errs by 10–25 % in the cooling time, as the
    red-sideband coupling saturates and reverses at ``n ≈ 3.67/η^2``.

`heating_rate` (quanta per time; one rate, or one per mode) adds the bath
operators ``\\sqrt{\\dot n} a_m`` and ``\\sqrt{\\dot n} a_m^†`` to every mode
with a non-zero rate. Laser-dephasing operators act
on the internal states only. Requires a static internal model (no beat
harmonics). Everything is sparse, and stripped to inverse `time_unit`.

!!! note
    Defined in the `LevelsQuantumToolboxExt` package extension, i.e. only once
    QuantumToolbox.jl is loaded.
"""
function motional_model end

"""
    mean_phonon_number(ρ, mm::MotionalModel; method = :direct)

Returns the mean phonon number of each mode from a density matrix on the space
of the [`MotionalModel`](@ref) (a vector for several modes, a number for one),
from the populations ``p_n`` of the reduced motional state, with `method`

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
    cooling_time(mm::MotionalModel; mode = nothing, eigvals = 12, dense_limit = 2500, weight_threshold = 0.1)
    cooling_time(H, c_ops; time_unit = u"µs", mode_sizes = nothing, kwargs...)

Returns the characteristic time of the slowest relaxation of the mean phonon
number (of the given `mode`, or of the total phonon number of all modes) under
the master equation of a [`MotionalModel`](@ref) — or of QuantumToolbox
operators `H`, `c_ops` on an internal ⊗ motional space, their motional factor
sizes inferred from the dimensions unless given — as the negative inverse real
part of the slowest Liouvillian eigenvalue, after the stationary state, whose
eigenmode has a weight ``|\\mathrm{Tr}[\\hat N v]| / ‖v‖`` in the phonon
number of at least `weight_threshold` times the largest such weight. Coherence
modes — which in the rate-equation limit relax at *half* the rate of ``\\bar n``
and would otherwise be picked up as the slowest, carrying only a weight of
order ``η²`` through the sideband coupling — are thereby skipped, so the
result matches the `τ_c` of [`cooling_rates`](@ref) where the adiabatic
elimination holds (and the Fock truncation is loose enough not to shape the
slowest mode). Liouvillians of dimension up to `dense_limit` are diagonalised
densely, larger ones with a shift-inverted Arnoldi iteration for `eigvals`
eigenpairs around zero. Where the relaxation is not a single exponential — a
hot ion whose Fock tail cools slowly — the asymptotic rate found here is not
the cooling time; [`IntegratedTransientSolver`](@ref) then gives the
state-dependent integral relaxation time instead.

!!! note
    Defined in the `LevelsQuantumToolboxExt` package extension, i.e. only once
    QuantumToolbox.jl is loaded.
"""
function cooling_time end

export RecoilLabel, HeatingLabel, MotionalModel
export populations, cooling_rates, motional_model, mean_phonon_number, cooling_time
