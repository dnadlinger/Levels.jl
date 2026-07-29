# The driven-transition model: rotating-frame Hamiltonian assembly from static
# field, harmonic drive terms, and the (possibly phase-modulated) laser coupling.

"""
A Hermitian drive term ``X \\cos(Ω t + φ)`` at the fundamental drive frequency
``Ω`` of a [`DrivenTransition`](@ref).
"""
struct HarmonicDrive{M<:AbstractMatrix}
    "Amplitude matrix ``X`` (Hermitian, in angular frequency units)."
    amplitude::M

    "Temporal phase ``φ``."
    phase::Float64
end

"""
    zeeman_drive(species, basis, B; phase = 0.0)

Returns the [`HarmonicDrive`](@ref) describing an ac magnetic field
``\\vec{B} \\cos(Ω t + φ)`` with the Cartesian amplitude vector `B`, via
[`Levels.zeeman_hamiltonian`](@ref).
"""
zeeman_drive(species, basis, B; phase=0.0) =
    HarmonicDrive(zeeman_hamiltonian(species, basis, B), Float64(phase))

"""
Sinusoidal phase modulation ``φ(t) = β_\\cos \\cos(Ω t) + β_\\sin \\sin(Ω t)`` of
the laser coupling of a [`DrivenTransition`](@ref) at its fundamental drive
frequency.

For excess micromotion at the trap rf frequency, ``β_\\cos`` is the in-phase and
``β_\\sin`` the quadrature modulation index.
"""
Base.@kwdef struct HarmonicPhaseModulation
    β_cos::Float64 = 0.0
    β_sin::Float64 = 0.0
end

"""
Returns a plain function mapping the time in seconds (`Float64`) to the laser
phase modulation in radians, for use inside the numeric engines.

`modulation` is either a [`HarmonicPhaseModulation`](@ref), or any callable taking
the (unitful) time to the phase; `Ω` is the fundamental drive frequency of the
[`DrivenTransition`](@ref) in angular units.
"""
function phase_function(modulation::HarmonicPhaseModulation, Ω)
    Ω_si = ustrip(u"s^-1", Ω)
    t -> modulation.β_cos * cos(Ω_si * t) + modulation.β_sin * sin(Ω_si * t)
end

phase_function(modulation, _) = t -> float(modulation(t * u"s"))

"""
Rotating-frame (optical RWA) model of a laser-probed transition between the
Zeeman states of two levels under periodic driving:

``H(t) = \\mathrm{diag}(\\texttt{frame}) - δ P_\\mathrm{upper}
+ \\sum_i X_i \\cos(Ω t + φ_i) + f(t) L + f(t)^* L^†``,

with ``f(t) = e^{-i φ_\\mathrm{mod}(t)} / 2`` the laser phase-modulation factor.
The laser detuning ``δ`` from the probed transition and the phase modulation
``φ_\\mathrm{mod}`` are supplied at evaluation time (cf.
[`full_hamiltonian`](@ref), [`sideband_amplitude`](@ref)), so one model instance
can serve e.g. a whole modulation-index scan.

Constructed via the `DrivenTransition(species, basis, probe; ...)` method; all
fields are in angular frequency units.
"""
Base.@kwdef struct DrivenTransition{
    B<:StateBasis,
    F<:Quantity,
    E<:Quantity,
    D<:HarmonicDrive,
    C<:Quantity,
}
    "The state basis the matrices are expressed in."
    basis::B

    "Fundamental angular frequency ``Ω`` of the periodic drive."
    drive_frequency::F

    "Rotating-frame energies at ``δ = 0`` (zero for the probed pair)."
    frame::Vector{E}

    "Harmonic drive terms ``X_i \\cos(Ω t + φ_i)``, block-diagonal in the levels."
    drives::Vector{D}

    "Laser coupling ``L`` in the upper⟨row|lower⟩⟨col| block; entries are carrier Rabi frequencies."
    coupling::Matrix{C}

    "Basis index of the probed lower state."
    lower::Int

    "Basis index of the probed upper state."
    upper::Int

    "Basis index range of the lower-level manifold."
    lower_range::UnitRange{Int}

    "Basis index range of the upper-level manifold."
    upper_range::UnitRange{Int}
end

"""
    DrivenTransition(species, basis, probe::Pair; drive_frequency, static_field,
                     coupling, drives = HarmonicDrive[])

Assembles the rotating-frame model for probing the `lower => upper` pair of
Zeeman states `probe`.

The rotating-frame energies are the Zeeman shifts in `static_field` (the static
magnetic flux density along the quantisation axis ẑ) relative to the probed pair.
`coupling` is the complex matrix of laser carrier Rabi frequencies (angular
units) in the upper⟨row|lower⟩⟨col| block, e.g. from
[`Levels.rabi_normalised`](@ref) for a beam geometry, or built directly for
couplings only present via mediator levels. `drives` are the harmonic terms at
`drive_frequency`, e.g. from [`zeeman_drive`](@ref).

The basis must consist of exactly the two probed levels; the harmonic drives
must be block-diagonal in them.
"""
function DrivenTransition(
    species::NoHyperfineOneElectronSpecies,
    basis::StateBasis,
    probe::Pair;
    drive_frequency,
    static_field,
    coupling,
    drives=HarmonicDrive[],
)
    lower_state = convert(StateSpec{NoHyperfineNumberSpec}, probe.first)
    upper_state = convert(StateSpec{NoHyperfineNumberSpec}, probe.second)
    if lower_state.level == upper_state.level
        throw(ArgumentError("Probed states must belong to two different levels"))
    end
    if Set(basis.levels) != Set([lower_state.level, upper_state.level])
        throw(
            ArgumentError(
                "Basis must consist of exactly the two probed levels " *
                "(spectator manifolds are not supported)",
            ),
        )
    end
    lower_range = staterange(basis, lower_state.level)
    upper_range = staterange(basis, upper_state.level)

    n = length(basis)
    if size(coupling) != (n, n)
        throw(ArgumentError("Laser coupling matrix does not match the basis size"))
    end
    if !all(
        iszero(coupling[i, k]) for i in 1:n for
        k in 1:n if !(i in upper_range && k in lower_range)
    )
        throw(
            ArgumentError(
                "Laser coupling must live in the upper⟨row|lower⟩⟨col| block only",
            ),
        )
    end

    for drive in drives
        if size(drive.amplitude) != (n, n)
            throw(ArgumentError("Drive amplitude matrix does not match the basis size"))
        end
        if !(
            all(iszero, drive.amplitude[lower_range, upper_range]) &&
            all(iszero, drive.amplitude[upper_range, lower_range])
        )
            throw(ArgumentError("Harmonic drives must be block-diagonal in the levels"))
        end
    end

    # Rotating-frame diagonal at δ = 0: Zeeman shifts relative to the probed pair.
    frame = map(collect(basis)) do state
        reference = state.level == lower_state.level ? lower_state : upper_state
        uconvert(
            u"µs^-1",
            zeeman_shift(species, state, static_field) -
            zeeman_shift(species, reference, static_field),
        )
    end

    # Store all matrices in one common unit; the concrete element type also keeps
    # the drive vector typed when it is empty.
    coupling = uconvert.(u"µs^-1", complex.(coupling))
    drives = HarmonicDrive{typeof(coupling)}[
        HarmonicDrive(uconvert.(u"µs^-1", complex.(d.amplitude)), d.phase) for
        d in drives
    ]

    DrivenTransition(;
        basis,
        drive_frequency,
        frame,
        drives,
        coupling,
        lower=stateindex(basis, lower_state),
        upper=stateindex(basis, upper_state),
        lower_range,
        upper_range,
    )
end

"""
    full_hamiltonian(dt::DrivenTransition, t; δ = 0, modulation = none)

Returns the full rotating-frame Hamiltonian matrix ``H(t)`` (in angular frequency
units) at the given (unitful) time, for laser detuning `δ` from the probed
transition and the given laser phase modulation (a
[`HarmonicPhaseModulation`](@ref), or any callable mapping the unitful time to
the phase in radians).
"""
function full_hamiltonian(
    dt::DrivenTransition,
    t;
    δ=zero(dt.drive_frequency),
    modulation=HarmonicPhaseModulation(),
)
    f = 0.5 * cis(-phase_function(modulation, dt.drive_frequency)(ustrip(u"s", t)))
    θ = uconvert(NoUnits, dt.drive_frequency * t)

    H = f .* dt.coupling .+ conj(f) .* Matrix(dt.coupling')
    for drive in dt.drives
        H .+= cos(θ + drive.phase) .* drive.amplitude
    end
    H .+= Diagonal(dt.frame)
    for i in dt.upper_range
        H[i, i] -= δ
    end
    H
end

export HarmonicDrive, zeeman_drive, HarmonicPhaseModulation, DrivenTransition
public full_hamiltonian
