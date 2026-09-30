"""
Multi-level optical Bloch equations: the Lindblad master equation of a chosen
set of levels of a species under a set of laser beams, assembled from the atomic
data — and, as a separate layer on top, its coupling to the harmonic motion of
the trapped ion for laser-cooling calculations.

A [`LaserScheme`](@ref) specifies the problem: the [`StateBasis`](@ref) of the
levels kept, the [`LaserBeam`](@ref)s (each a [`RelativeFrequency`](@ref)
detuning from a zero-field line centre, an intensity, a polarisation, a
direction and optionally a linewidth) and the static magnetic field.
[`lindblad_model`](@ref) turns it into a [`LindbladModel`](@ref): the
rotating-frame Hamiltonian — each fine-structure level (or, on request, each
state) given one frame frequency (cf. [`RotatingFrame`](@ref)) such that every
beam coupling is time-independent where the beam graph allows it, with beat
notes kept as harmonic terms otherwise — plus the spontaneous-emission and
laser-phase-diffusion jump operators. Everything is species-kind agnostic: for a
[`Levels.HyperfineOneElectronSpecies`](@ref) the basis states denote the
adiabatically-labelled eigenstates at the static field (cf.
[`hyperfine_manifold`](@ref)), with exact at-field coupling and decay
amplitudes.

The motional layer never enters that model: a [`MotionalMode`](@ref) (frequency
and direction) combines with a `LindbladModel` into a
[`MotionalCoupling`](@ref) — the projected Lamb–Dicke factors per beam, the
first-order sideband Hamiltonian and the recoil jump operators, one instance
per mode — consumed by the QuantumToolbox-based solvers of the
`LevelsQuantumToolboxExt` package extension (steady states, adiabatic-elimination
cooling rates, full internal ⊗ motional models).

All quantities are unitful and, as in [`Levels.PeriodicDynamics`](@ref),
normalised to µs⁻¹ (jump operators to µs⁻¹ᐟ²).
"""
module OpticalBloch

using LinearAlgebra
using StaticArrays: SVector
using Unitful

using ..Levels:
    Levels,
    HyperfineNumberSpec,
    LevelSpec,
    NoHyperfineNumberSpec,
    OneElectronSpecies,
    RelativeFrequency,
    StateBasis,
    StateSpec,
    amplitude_evaluator,
    coupling_matrix,
    einstein_a,
    fine_structure,
    hyperfine_manifold,
    hyperfine_shift,
    multipole_rank,
    parse_level,
    photon_energy,
    state_energy,
    stateindex,
    transition_frequency,
    zeeman_shift

const ANGULAR_UNIT = u"µs^-1"
const JUMP_UNIT = u"µs^(-1/2)"

include("scheme.jl")
include("frame.jl")
include("lindblad.jl")
include("motion.jl")
include("quantum_toolbox.jl")

end # module
