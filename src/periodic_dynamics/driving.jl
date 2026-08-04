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
    zeeman_drives(species, basis, b; phase = 0.0)

Returns the [`HarmonicDrive`](@ref)s describing an ac magnetic field given as a
complex Cartesian amplitude phasor `b`,
``\\vec{B}(t) = \\mathrm{Re}[\\vec{b} \\, e^{-i (Ω t + φ)}]``: one drive
``\\mathrm{Re}(\\vec{b}) \\cos(Ω t + φ)`` plus, for genuinely complex `b`, a
second ``\\mathrm{Im}(\\vec{b}) \\cos(Ω t + φ - π/2)`` (identically zero
components are dropped). A real `b` is a linearly polarised field; e.g.
``\\vec{b} = B_0/\\sqrt{2} \\, (1, -i, 0)`` rotates in the x–y plane.
"""
function zeeman_drives(species, basis, b; phase=0.0)
    drives = HarmonicDrive[]
    re = real.(collect(b))
    if !all(iszero, re)
        push!(drives, zeeman_drive(species, basis, re; phase=Float64(phase)))
    end
    im_part = imag.(collect(b))
    if !all(iszero, im_part)
        push!(
            drives,
            zeeman_drive(species, basis, im_part; phase=Float64(phase) - π / 2),
        )
    end
    drives
end

"""
Sinusoidal phase modulation
``φ(t) = β_\\mathrm{cos} \\cos(Ω t) + β_\\mathrm{sin} \\sin(Ω t)`` of the laser coupling
of a [`DrivenTransition`](@ref) at its fundamental drive frequency.

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

```math
H(t) = \\mathrm{diag}(\\texttt{frame}) - δ P_\\mathrm{upper}
+ \\sum_i X_i \\cos(Ω t + φ_i) + f(t) L + f(t)^* L^†,
```

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

# Validates the coupling and drive matrices of a DrivenTransition against the
# basis size and the two manifold index ranges.
function validate_driven_blocks(n, coupling, drives, lower_range, upper_range)
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
            throw(
                ArgumentError(
                    "Harmonic drives must be block-diagonal in the manifolds",
                ),
            )
        end
    end
end

# Normalises all matrices of a DrivenTransition to one common unit; the concrete
# element type also keeps the drive vector typed when it is empty.
function normalise_driven_matrices(coupling, drives)
    coupling = uconvert.(u"µs^-1", complex.(coupling))
    drives = HarmonicDrive{typeof(coupling)}[
        HarmonicDrive(uconvert.(u"µs^-1", complex.(d.amplitude)), d.phase) for
        d in drives
    ]
    coupling, drives
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

For a [`HyperfineOneElectronSpecies`](@ref), the basis must consist of the two
**complete** fine-structure manifolds the probed hyperfine states belong to
(all ``F``, all ``m_F``, contiguous per manifold — e.g. from
`StateBasis(species, "S_1/2", "D_5/2")`). The basis states then denote the
**adiabatically-labelled eigenstates** of the hyperfine + Zeeman Hamiltonian at
`static_field` (cf. [`Levels.hyperfine_manifold`](@ref)): the rotating-frame
energies are the exact manifold eigen-energies, and `coupling` and `drives` —
supplied in the coupled zero-field ``|F, m_F⟩`` basis, e.g. from
[`Levels.quadrupole_couplings`](@ref) and [`zeeman_drives`](@ref) — are rotated
into that eigenbasis internally, which accounts exactly for the ``F`` mixing at
the working field.
"""
function DrivenTransition(
    species::NoHyperfineOneElectronSpecies,
    basis::StateBasis{NoHyperfineNumberSpec},
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

    validate_driven_blocks(length(basis), coupling, drives, lower_range, upper_range)

    # Rotating-frame diagonal at δ = 0: Zeeman shifts relative to the probed pair.
    frame = map(collect(basis)) do state
        reference = state.level == lower_state.level ? lower_state : upper_state
        uconvert(
            u"µs^-1",
            zeeman_shift(species, state, static_field) -
            zeeman_shift(species, reference, static_field),
        )
    end

    coupling, drives = normalise_driven_matrices(coupling, drives)

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

# Canonicalises a probed-state specification and checks it names a hyperfine
# (F) level.
function hyperfine_probe_state(state::StateSpec)
    spec = parse_level(state.level)
    if !(spec isa HyperfineNumberSpec)
        throw(
            ArgumentError(
                "Probed state '$state' must specify a hyperfine (F) level " *
                "for a hyperfine species",
            ),
        )
    end
    StateSpec(spec, state.m)
end

function DrivenTransition(
    species::HyperfineOneElectronSpecies,
    basis::StateBasis{HyperfineNumberSpec},
    probe::Pair;
    drive_frequency,
    static_field,
    coupling,
    drives=HarmonicDrive[],
)
    lower_state = hyperfine_probe_state(probe.first)
    upper_state = hyperfine_probe_state(probe.second)
    fs_lower = fine_structure(lower_state.level)
    fs_upper = fine_structure(upper_state.level)
    if fs_lower == fs_upper
        throw(
            ArgumentError(
                "Probed states must belong to two different fine-structure manifolds",
            ),
        )
    end
    lower_range = staterange(basis, fs_lower)
    upper_range = staterange(basis, fs_upper)
    n = length(basis)
    if length(lower_range) + length(upper_range) != n
        throw(
            ArgumentError(
                "Basis must consist of exactly the two probed manifolds " *
                "(spectator manifolds are not supported)",
            ),
        )
    end
    for (fs, range) in ((fs_lower, lower_range), (fs_upper, upper_range))
        if length(range) != Int(2 * species.nuclear_spin + 1) * Int(2 * fs.j + 1)
            throw(
                ArgumentError(
                    "Basis must contain the complete '$fs' manifold — the " *
                    "eigenbasis rotation needs the full state space",
                ),
            )
        end
    end

    validate_driven_blocks(n, coupling, drives, lower_range, upper_range)

    m_lower = hyperfine_manifold(species, fs_lower, static_field)
    m_upper = hyperfine_manifold(species, fs_upper, static_field)

    # Permutations from the (possibly reordered) user basis to the canonical
    # manifold bases; the eigen-solutions are label-aligned, so this reorders
    # both the coupled-basis components and the adiabatic labels consistently.
    perm_lower = [stateindex(m_lower.basis, basis[i]) for i in lower_range]
    perm_upper = [stateindex(m_upper.basis, basis[i]) for i in upper_range]
    v_lower = m_lower.states[perm_lower, perm_lower]
    v_upper = m_upper.states[perm_upper, perm_upper]

    # Rotating-frame diagonal at δ = 0: exact eigen-energies relative to the
    # probed pair.
    frame = Vector{eltype(m_lower.energies)}(undef, n)
    frame[lower_range] =
        m_lower.energies[perm_lower] .-
        m_lower.energies[stateindex(m_lower.basis, lower_state)]
    frame[upper_range] =
        m_upper.energies[perm_upper] .-
        m_upper.energies[stateindex(m_upper.basis, upper_state)]

    # Rotate the coupled-basis coupling and drives into the field eigenbasis.
    coupling = complex.(coupling)
    rotated_coupling = zero(coupling)
    rotated_coupling[upper_range, lower_range] =
        v_upper' * coupling[upper_range, lower_range] * v_lower
    rotated_drives = map(drives) do d
        amplitude = complex.(d.amplitude)
        rotated = zero(amplitude)
        rotated[lower_range, lower_range] =
            v_lower' * amplitude[lower_range, lower_range] * v_lower
        rotated[upper_range, upper_range] =
            v_upper' * amplitude[upper_range, upper_range] * v_upper
        HarmonicDrive(rotated, d.phase)
    end

    coupling, drives = normalise_driven_matrices(rotated_coupling, rotated_drives)

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

export HarmonicDrive,
    zeeman_drive, zeeman_drives, HarmonicPhaseModulation, DrivenTransition
public full_hamiltonian
