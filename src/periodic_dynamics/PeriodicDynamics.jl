"""
Periodic driving of laser-probed atomic transitions: ac Zeeman coupling,
micromotion phase modulation, sideband Rabi frequencies and resonance shifts.

An [`Eigenstructure`](@ref Levels.Eigenstructure) supplies the static atom.
[`DrivenTransition`](@ref) adds the probe laser and harmonic drives in the
optical rotating-wave approximation. Its matrices use the structure's field
eigenbasis, including field-induced hyperfine mixing.

[`sideband_rabi`](@ref) offers two methods: [`DressedStates`](@ref) dresses
the manifolds with the periodic drive and treats the laser to first order;
[`Monodromy`](@ref) propagates the full model and finds an avoided-crossing
gap. [`quasienergy_shift`](@ref) evaluates laser-free shifts with
[`PerturbationTheory`](@ref) or nonperturbative [`Floquet`](@ref) dressing.

Cf. Joshi et al., "Characterization of ion-trap-induced ac magnetic fields",
Phys. Rev. A **110**, 063101 (2024),
[doi:10.1103/PhysRevA.110.063101](https://doi.org/10.1103/PhysRevA.110.063101).
"""
module PeriodicDynamics

using LinearAlgebra
using Unitful

using ..Levels:
    Levels,
    Eigenstructure,
    StateBasis,
    StateSpec,
    fine_structure,
    iscomplete,
    offsets,
    parse_level,
    stateindex,
    staterange,
    validate_state_kind,
    zeeman_hamiltonian

include("methods.jl")
include("driving.jl")
include("floquet.jl")
include("quasienergy.jl")
include("tdse.jl")
include("quantum_toolbox.jl")

end # module
