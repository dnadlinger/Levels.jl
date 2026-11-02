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
The internal ⊗ Fock master equation built by [`motional_model`](@ref).

`H` and `c_ops` are QuantumToolbox operators with the internal factor first,
then one Fock space per mode (last mode fastest). `H` uses inverse `time_unit`
and jump operators inverse square-root `time_unit`; their entries are plain
numbers. `c_labels` identifies internal decay, dephasing, recoil and heating
operators. `num_states`, `num_fock`, `modes` and `lamb_dicke_order` record the
layout and displacement expansion.

With QuantumToolbox loaded, the model supports steady states, time evolution,
phonon populations and the relaxation methods of [`cooling_metrics`](@ref).
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
    motional_model(model::LindbladModel, mc::MotionalCoupling; num_fock, lamb_dicke_order = 1, heating_rate = 0, quadrature_order = nothing, time_unit = u"µs")
    motional_model(model::LindbladModel, mcs::AbstractVector{<:MotionalCoupling}; num_fock, …) -> MotionalModel

Builds an internal ⊗ Fock [`MotionalModel`](@ref), with `num_fock` states per
mode (one number or one per mode). For the full displacement operators,

```math
H = H_\\mathrm{frame}⊗1 + ∑_m ω_m a_m^†a_m
  + \\frac12∑_b(C_b⊗D_b+\\mathrm{h.c.}), \\qquad
D_b = e^{i∑_m η_{b,m}(a_m+a_m^†)}.
```

``H_\\mathrm{frame}`` is the diagonal part of `model.hamiltonian`, excluding
its beam couplings. ``C_b`` is each beam's Rabi-frequency matrix and
``η_{b,m}`` its projected Lamb–Dicke factor. Spontaneous emission contributes
``L⊗e^{-i\\vec k⋅\\vec r}`` averaged over its radiation pattern (Kulosa et al.,
New J. Phys. **25**, 053008 (2023), Eqs. (6), (7), (A.3)).

`lamb_dicke_order` selects the displacement approximation:

- `1` (default) keeps the linear sideband Hamiltonian. Each decay operator
  acts as ``L⊗1`` with separate recoil operators
  ``\\sqrt{α_m}η_{0,m}L⊗(a_m+a_m^†)``. This model reproduces the recoil heating
  ``α_mη_{0,m}^2`` per emission and needs only the second emission moment;
  it accepts a constant `recoil_moment` and drops cross-mode recoil correlations.
- A finite integer ``k ≥ 2`` truncates beam and recoil displacements at total
  order ``k``. Recoil uses [`emission_rule`](@ref) with default
  `quadrature_order = k + 1`.
- `Inf` uses the full displacement matrix elements of
  [`displacement_elements`](@ref). Each decay is replaced by quadrature
  operators ``\\sqrt{w_j}L⊗D(-η_0\\vec h_j)``. The default `quadrature_order`
  is 12. Several modes use joint projections of emission directions from a
  product rule on the sphere, retaining their recoil correlations.

Orders above 1 require emission distributions from `recoil_moment = :exact`
or `:isotropic`. For one mode, `quadrature_order` is the number of emission
nodes; for several modes it is the polar order of the spherical product rule.
Check convergence in Fock size and quadrature order; for finite expansions,
also check the Lamb–Dicke order.

`heating_rate` (one rate or one per mode) adds
``\\sqrt{\\dot n}a_m`` and ``\\sqrt{\\dot n}a_m^†``. Laser phase diffusion acts on
internal states only. Requires a static internal model, without beat harmonics.
Operators are sparse; Hamiltonians are stripped to inverse `time_unit` and
jumps to inverse square-root `time_unit`.

!!! note
    Defined when QuantumToolbox.jl is loaded.
"""
function motional_model end

"""
    mean_phonon_number(ρ, mm::MotionalModel; method = :direct)

Returns the mean phonon number of each mode (a number for one mode, a vector
for several) from its reduced populations ``p_n``, ``n=0,…,N-1``, where
``N`` is that mode's `num_fock`.

- `:direct` returns ``∑_n n p_n``. Increase the Fock truncation until this
  expectation converges.
- `:thermal_ratio` returns ``r/(1-r)`` with
  ``r = (∑_{n=1}^{N-1}p_n)/(∑_{n=0}^{N-2}p_n)``. For a truncated thermal
  distribution ``p_n∝r^n`` this recovers the untruncated mean at any ``N≥2``.
  For a nonthermal state it is an estimator: as the truncation grows it tends
  to ``1/p_0-1``, which need not equal the mean.

!!! note
    Defined when QuantumToolbox.jl is loaded.
"""
function mean_phonon_number end

# --- Common interface of the cooling solvers ----------------------------------------

"""
A method for estimating the laser-cooling performance of a mode, dispatched on
by [`cooling_metrics`](@ref): [`AdiabaticElimination`](@ref),
[`IntegratedTransient`](@ref) or [`LiouvillianSpectrum`](@ref).
"""
abstract type CoolingMethod end

"""
    AdiabaticElimination()

Selects the rate equations obtained by adiabatically eliminating the internal
dynamics: steady-state occupation and an exponential time constant from rates
of order ``η^2``. Useful for parameter searches when the approximations hold;
full-model methods resolve saturation and nonthermal motion beyond that regime.
Requires the Lamb–Dicke regime and weak sideband coupling, with internal
relaxation fast compared with the induced motional rates. The model must
have no beat-note harmonics.

```math
A_± = 2 \\mathrm{Re} S(±ω_m) + D, \\qquad
S(ω) = \\mathrm{Tr}\\left[H_\\mathrm{sb}(-\\mathcal L+iω)^{-1}
    (H_\\mathrm{sb}ρ_∞)\\right]
= \\int_0^∞ e^{-iωτ}⟨H_\\mathrm{sb}(τ)H_\\mathrm{sb}(0)⟩\\,dτ,
```

where ``\\mathcal L`` is the internal Liouvillian and ``ρ_∞`` its steady state.
Spontaneous-emission recoil adds
``D = ∑_k α_k η_{0,k}^2 \\mathrm{Tr}(L_k^†L_kρ_∞)``; see
[`MotionalCoupling`](@ref), Cirac et al., Phys. Rev. A **46**, 2668 (1992),
and Morigi, Phys. Rev. A **67**, 033402 (2003).

With this Fourier sign, cooling samples ``S(-ω_m)``. In the weak-excitation
two-level limit,
``A_∓-D = η^2 Γ|Ω|^2/[Γ^2+4(δ±ω_m)^2]`` with laser-minus-atomic detuning
``δ``: a red sideband at ``δ=-ω_m`` enhances ``A_-``.

``A_+`` raises and ``A_-`` lowers the phonon number:

```math
\\dot{\\bar n} = A_+(\\bar n+1)-A_-\\bar n+\\dot n_\\mathrm{heat}, \\qquad
\\bar n_∞ = \\frac{A_++\\dot n_\\mathrm{heat}}{A_--A_+}, \\qquad
τ_c = \\frac{1}{A_--A_+}.
```

The [`CoolingMetrics`](@ref) fields `nbar` and `τ_c` are `NaN` when
``A_- ≤ A_+``. Several modes share the internal steady state and factorisation
and are independent at this order.
"""
struct AdiabaticElimination <: CoolingMethod end

"""
    IntegratedTransient(; initial_nbar, num_fock = fock_truncation.(initial_nbar), lamb_dicke_order = Inf, bandwidth = 4, quadrature_order = nothing)

Selects the integrated-transient solution of a [`motional_model`](@ref),
restricted to motional coherences within `bandwidth` by
[`BandLiouvillian`](@ref). Returns the steady occupation of this truncated
model and the integral relaxation time from the internal steady state ⊗ a
thermal state of `initial_nbar` (one value or one per mode).

`num_fock` defaults to a thermal-tail tolerance of ``10^{-4}`` via
[`fock_truncation`](@ref). `bandwidth = nothing` retains all coherences;
`lamb_dicke_order` and `quadrature_order` pass to `motional_model`.
The result carries an [`IntegratedTransientSolver`](@ref) for relaxation
times from other initial states without repeating the factorisation.
"""
struct IntegratedTransient <: CoolingMethod
    initial_nbar::Union{Float64,Vector{Float64}}
    num_fock::Union{Nothing,Int,Vector{Int}}
    lamb_dicke_order::Float64
    bandwidth::Union{Nothing,Int}
    quadrature_order::Union{Nothing,Int}
end

function IntegratedTransient(;
    initial_nbar,
    num_fock=nothing,
    lamb_dicke_order::Real=Inf,
    bandwidth=4,
    quadrature_order=nothing,
)
    validate_lamb_dicke_order(lamb_dicke_order)
    nb = initial_nbar isa Real ? Float64(initial_nbar) : collect(Float64, initial_nbar)
    all(>=(0), nb) || throw(ArgumentError("initial_nbar must be non-negative"))
    nf =
        isnothing(num_fock) ? nothing :
        num_fock isa Integer ? Int(num_fock) : collect(Int, num_fock)
    IntegratedTransient(
        nb,
        nf,
        Float64(lamb_dicke_order),
        isnothing(bandwidth) ? nothing : Int(bandwidth),
        isnothing(quadrature_order) ? nothing : Int(quadrature_order),
    )
end

"""
    LiouvillianSpectrum(; num_fock, lamb_dicke_order = 1, num_eigenpairs = 12, dense_limit = 2500, weight_threshold = 0.1)

Selects the full internal ⊗ Fock master equation ([`motional_model`](@ref) with
`num_fock` states per mode, to the given `lamb_dicke_order`), solved for its
stationary state and slowest eigenmodes from one trace-row factorisation
([`relaxation_modes`](@ref)). The occupation comes from the former. The time
constant is the asymptotic phonon-number relaxation time
``-1/\\mathrm{Re}\\,λ`` of the slowest among the `num_eigenpairs` slowest modes
whose phonon-number weight ``|\\mathrm{Tr}(\\hat N_m v)|/‖v‖`` is at least
`weight_threshold` times the largest weight among them; it is state-independent
like the adiabatic-elimination one, but without its approximation.

The threshold passes over coherence modes, which carry a phonon-number weight
of order ``η^2`` only (the slowest of them decays at half the cooling rate).
Increase `num_eigenpairs` and vary the threshold to check the selection; an
error is raised if no mode qualifies. `dense_limit` is the Liouvillian dimension
up to which dense diagonalisation is used. A negative time signals a growing
mode of the truncated model.

For a transient with several timescales this asymptotic time differs from the
initial-state-dependent time of [`IntegratedTransient`](@ref).
"""
struct LiouvillianSpectrum <: CoolingMethod
    num_fock::Union{Int,Vector{Int}}
    lamb_dicke_order::Float64
    bandwidth::Union{Nothing,Int}
    num_eigenpairs::Int
    dense_limit::Int
    weight_threshold::Float64
end

function LiouvillianSpectrum(;
    num_fock,
    lamb_dicke_order::Real=1,
    bandwidth=nothing,
    num_eigenpairs::Integer=12,
    dense_limit::Integer=2500,
    weight_threshold::Real=0.1,
)
    validate_lamb_dicke_order(lamb_dicke_order)
    LiouvillianSpectrum(
        num_fock isa Integer ? Int(num_fock) : collect(Int, num_fock),
        Float64(lamb_dicke_order),
        isnothing(bandwidth) ? nothing : Int(bandwidth),
        num_eigenpairs,
        dense_limit,
        weight_threshold,
    )
end

"""
The cooling performance of one mode as estimated by [`cooling_metrics`](@ref)
with the given [`CoolingMethod`](@ref):

  - `nbar`: the steady-state mean phonon number (`NaN` where the
    adiabatic-elimination rates predict heating);
  - `τ_c`: the cooling time constant (unitful) — the exponential time constant
    ``1/(A_- - A_+)`` for `AdiabaticElimination` (`NaN` if heated), the
    integral relaxation time of the phonon number from the method's initial
    state for `IntegratedTransient`, the slowest phonon-number relaxation for
    `LiouvillianSpectrum`;
  - `A_plus`, `A_minus`: the rate coefficients (`AdiabaticElimination` only,
    `nothing` otherwise);
  - `model`: the [`MotionalModel`](@ref) built for the estimate (full-model
    methods only, `nothing` otherwise);
  - `solver`: the [`IntegratedTransientSolver`](@ref) (`IntegratedTransient`
    only, `nothing` otherwise), for relaxation times from other initial
    states.

For several modes, `cooling_metrics` returns one `CoolingMetrics` per mode;
the full-model methods then share one joint model and solver.
"""
struct CoolingMetrics{M<:CoolingMethod,T,R,F,S}
    method::M
    nbar::Float64
    τ_c::T
    A_plus::R
    A_minus::R
    model::F
    solver::S
end

function Base.show(io::IO, m::CoolingMetrics)
    print(
        io,
        "CoolingMetrics(",
        nameof(typeof(m.method)),
        ": n̄ = ",
        round(m.nbar; sigdigits=4),
        ", τ_c = ",
        m.τ_c,
        ")",
    )
end

"""
    cooling_metrics(model::LindbladModel, mc::MotionalCoupling, method::CoolingMethod; heating_rate = 0, ρ = steadystate(model), time_unit = u"µs") -> CoolingMetrics
    cooling_metrics(model::LindbladModel, mcs::AbstractVector{<:MotionalCoupling}, method; …) -> Vector{CoolingMetrics}

Returns the steady-state occupation and the cooling time constant of the
mode(s) under the laser scheme of `model`, as a [`CoolingMetrics`](@ref) per
mode, by the chosen [`CoolingMethod`](@ref): [`AdiabaticElimination`](@ref)
(rate equations on the internal steady state), [`IntegratedTransient`](@ref)
(full model, integral relaxation time from a thermal initial state) or
[`LiouvillianSpectrum`](@ref) (full model, asymptotic relaxation rate).
`heating_rate` (quanta per time) enters the rate equations as
``\\dot n_\\mathrm{heat}`` and the full models as a thermal bath (cf.
[`motional_model`](@ref)); `ρ` is the internal steady state, reused by the rate
equations and the internal factor of the full models' initial state. Rates and
times are in `time_unit`.

!!! note
    Defined in the `LevelsQuantumToolboxExt` package extension, i.e. only once
    QuantumToolbox.jl is loaded.
"""
function cooling_metrics end

Base.show(io::IO, model::MotionalModel) = print(
    io,
    "MotionalModel(",
    model.num_states,
    " internal states, num_fock=",
    model.num_fock,
    ", time_unit=",
    model.time_unit,
    ")",
)

export RecoilLabel, HeatingLabel, MotionalModel
export populations, motional_model, mean_phonon_number
export CoolingMethod, AdiabaticElimination, IntegratedTransient, LiouvillianSpectrum
export CoolingMetrics, cooling_metrics
