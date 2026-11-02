"""
Multi-level Lindblad master equations and their coupling to trapped-ion motion
for laser-cooling calculations.

[`lindblad_model`](@ref) combines an
[`Eigenstructure`](@ref Levels.Eigenstructure) with a [`LaserScheme`](@ref):
the structure supplies the field eigenstates, the scheme supplies beams and
a rotating-frame choice. The resulting [`LindbladModel`](@ref) contains the
Hamiltonian, spontaneous-emission operators and laser phase diffusion. Beam
couplings that cannot be made static in the chosen frame remain as beat-note
harmonics.

[`motional_coupling`](@ref) adds a [`MotionalMode`](@ref), giving projected
Lamb–Dicke factors, sideband couplings and recoil data. With QuantumToolbox
loaded, [`motional_model`](@ref) builds the internal ⊗ Fock master equation
and [`cooling_metrics`](@ref) compares rate-equation, integrated-transient and
Liouvillian-spectrum estimates. [`PhononRateEquation`](@ref) keeps the rate
description valid for saturated sidebands, without a joint model.

The solver-independent [`BandLiouvillian`](@ref) restricts motional
coherences for numerical calculations. [`IntegratedTransientSolver`](@ref)
uses a sparse linear solve for steady states and integral relaxation times.
Hamiltonians use angular-frequency units and jump operators square-root rate
units; see [Conventions and examples](@ref conventions).
"""
module OpticalBloch

using ArnoldiMethod: partialeigen, partialschur
using LinearAlgebra
using FastGaussQuadrature: gausslegendre
using Interpolations: linear_interpolation
using OrdinaryDiffEqVerner: Vern7
using SciMLBase: ContinuousCallback, ODEProblem, solve, successful_retcode, terminate!
using SparseArrays:
    SparseArrays,
    SparseMatrixCSC,
    dropzeros!,
    droptol!,
    nonzeros,
    nzrange,
    rowvals,
    sparse,
    spdiagm,
    spzeros
using StaticArrays: SVector
using Unitful

using ..Levels:
    Levels,
    Eigenstructure,
    HyperfineNumberSpec,
    LevelSpec,
    NoHyperfineNumberSpec,
    OneElectronSpecies,
    RelativeFrequency,
    StateBasis,
    StateSpec,
    amplitude_lookup,
    coupling_matrix,
    einstein_a,
    lifetime,
    fine_structure,
    hyperfine_shift,
    multipole_rank,
    offsets,
    parse_level,
    photon_energy,
    stateindex,
    transition_frequency

const ANGULAR_UNIT = u"rad/µs"
const JUMP_UNIT = u"µs^(-1/2)"

include("scheme.jl")
include("frame.jl")
include("lindblad.jl")
include("motion.jl")
include("liouvillian.jl")
include("quantum_toolbox.jl")
include("rate_equation.jl")

end # module
