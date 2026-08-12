"""
Simulation of optical transitions between Zeeman manifolds under periodic driving,
e.g. by the rf fields of a Paul trap: ac Zeeman coupling, laser phase modulation
from excess micromotion, and the resulting sidebands, shifts, and modifications to
the effective Rabi frequencies.

The model is assembled in the frame rotating at the laser frequency (optical RWA):
a [`DrivenTransition`](@ref) collects the rotating-frame energies for a static
magnetic field, any number of harmonic drive terms at the fundamental drive
frequency (e.g. an ac magnetic field via [`zeeman_drive`](@ref)), and the laser
coupling matrix, optionally phase-modulated (e.g. [`HarmonicPhaseModulation`](@ref)
for excess micromotion). Two engines evaluate it: fast dressed-state Floquet
perturbation theory ([`dress`](@ref), [`sideband_amplitude`](@ref)), and exact
monodromy-matrix propagation of the time-dependent Schrödinger equation
([`exact_sideband`](@ref), [`stroboscopic_populations`](@ref)) as a
nonperturbative cross-check.

The ac Zeeman shifts themselves — with or without a probe laser — are available
directly through the second-order perturbation sum [`ac_zeeman_shift`](@ref)
(with [`ac_zeeman_contributions`](@ref) diagnostics and near-resonance warnings)
and its nonperturbative counterpart [`floquet_zeeman_shift`](@ref); for
hyperfine species these operate on the exact
[`Levels.hyperfine_manifold`](@ref) eigenstates, whose ``F`` mixing the
hyperfine [`DrivenTransition`](@ref) method likewise accounts for by rotating
couplings and drives into the field eigenbasis.

Cf. Joshi et al., "Characterization of ion-trap-induced ac magnetic fields",
Phys. Rev. A **110**, 063101 (2024),
[doi:10.1103/PhysRevA.110.063101](https://doi.org/10.1103/PhysRevA.110.063101).
"""
module PeriodicDynamics

using LinearAlgebra
using Unitful

using ..Levels:
    Levels,
    HyperfineManifold,
    HyperfineNumberSpec,
    HyperfineOneElectronSpecies,
    NoHyperfineNumberSpec,
    NoHyperfineOneElectronSpecies,
    StateBasis,
    StateSpec,
    eigenbasis_transform,
    fine_structure,
    hyperfine_manifold,
    parse_level,
    state_energy,
    stateindex,
    staterange,
    zeeman_hamiltonian,
    zeeman_shift

include("driving.jl")
include("floquet.jl")
include("ac_zeeman.jl")
include("tdse.jl")
include("quantum_toolbox.jl")

end # module
