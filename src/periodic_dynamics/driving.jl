# The driven-transition model: rotating-frame Hamiltonian assembly from the
# structure's eigen-energies, harmonic drive terms, and the (possibly
# phase-modulated) laser coupling.

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
Returns the ``e^{+iΩt}`` Fourier component ``H_1 = Σ_i ½ e^{iφ_i} X_i``
of the given [`HarmonicDrive`](@ref)s ``Σ_i X_i \\cos(Ω t + φ_i)``, restricted
to the basis index range `r`, so that the drive term reads
``H_1 e^{iΩt} + H_1^† e^{-iΩt}`` (the convention of [`dress_manifold`](@ref)
and of the perturbative [`quasienergy_shift`](@ref)).
"""
function harmonic_block(drives, r)
    H1 = zeros(ComplexF64, length(r), length(r)) .* u"rad/µs"
    for drive in drives
        H1 .+= (0.5 * cis(drive.phase)) .* drive.amplitude[r, r]
    end
    H1
end

"""
    zeeman_drive(s::Eigenstructure, b; phase = 0.0)

Returns the [`HarmonicDrive`](@ref) describing an ac magnetic field
``\\vec{b} \\cos(Ω t + φ)`` with the real Cartesian amplitude vector `b`, via
[`Levels.zeeman_hamiltonian`](@ref)`(s, b)`: the drive matrix is expressed in
the field eigenbasis of the structure, indexed by its basis, like every
matrix of a [`DrivenTransition`](@ref) built on `s`.
"""
zeeman_drive(s::Eigenstructure, b; phase=0.0) =
    HarmonicDrive(zeeman_hamiltonian(s, b), Float64(phase))

"""
    zeeman_drives(s::Eigenstructure, b; phase = 0.0)

Returns the [`HarmonicDrive`](@ref)s describing an ac magnetic field given as a
complex Cartesian amplitude phasor `b`,
``\\vec{B}(t) = \\mathrm{Re}[\\vec{b} \\, e^{-i (Ω t + φ)}]``: one drive
``\\mathrm{Re}(\\vec{b}) \\cos(Ω t + φ)`` plus, for genuinely complex `b`, a
second ``\\mathrm{Im}(\\vec{b}) \\cos(Ω t + φ - π/2)`` (identically zero
components are dropped). A real `b` is a linearly polarised field; e.g.
``\\vec{b} = B_0/\\sqrt{2} \\, (1, -i, 0)`` rotates in the x–y plane. The
matrices are in the field eigenbasis of `s` (cf. [`zeeman_drive`](@ref)).
"""
function zeeman_drives(s::Eigenstructure, b; phase=0.0)
    drives = HarmonicDrive[]
    re = real.(collect(b))
    if !all(iszero, re)
        push!(drives, zeeman_drive(s, re; phase=Float64(phase)))
    end
    im_part = imag.(collect(b))
    if !all(iszero, im_part)
        push!(drives, zeeman_drive(s, im_part; phase=Float64(phase) - π / 2))
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
    Ω_si = ustrip(u"rad/s", Ω)
    t -> modulation.β_cos * cos(Ω_si * t) + modulation.β_sin * sin(Ω_si * t)
end

phase_function(modulation, _) = t -> float(modulation(t * u"s"))

"""
Rotating-frame model of a laser-probed transition between two manifolds, in the
optical rotating-wave approximation:

```math
H(t) = \\mathrm{diag}(\\texttt{frame}) - δ P_\\mathrm{upper}
+ \\sum_i X_i \\cos(Ω t + φ_i) + f(t) L + f(t)^* L^†,
\\qquad f(t) = e^{-iφ_\\mathrm{mod}(t)}/2.
```

`frame`, drive amplitudes `X_i`, laser coupling `L` and detuning `δ` are in
angular-frequency units. Detuning from the probed component and phase
modulation are supplied at evaluation time through [`full_hamiltonian`](@ref)
or [`sideband_rabi`](@ref), so one model can serve a parameter scan.

`DrivenTransition(s, probe; …)` constructs the model in the structure's field
eigenbasis. The keyword constructor
`DrivenTransition(; basis, drive_frequency, frame, drives, coupling, lower,
upper, lower_range, upper_range)` accepts manually assembled data without
validation; all arrays must use the declared basis and ordering.
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
function validate_driven_blocks(num_states, coupling, drives, lower_range, upper_range)
    if size(coupling) != (num_states, num_states)
        throw(ArgumentError("Laser coupling matrix does not match the basis size"))
    end
    if !all(
        iszero(coupling[i, k]) for i in 1:num_states for
        k in 1:num_states if !(i in upper_range && k in lower_range)
    )
        throw(
            ArgumentError(
                "Laser coupling must live in the upper⟨row|lower⟩⟨col| block only",
            ),
        )
    end

    for drive in drives
        if size(drive.amplitude) != (num_states, num_states)
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
    coupling = uconvert.(u"rad/µs", complex.(coupling))
    drives = HarmonicDrive{typeof(coupling)}[
        HarmonicDrive(uconvert.(u"rad/µs", complex.(d.amplitude)), d.phase) for
        d in drives
    ]
    coupling, drives
end

"""
    DrivenTransition(s::Eigenstructure, probe::Pair; drive_frequency, coupling,
                     drives = HarmonicDrive[])

Assembles a periodic model probing the `lower => upper` states of `s`.
The static frame diagonal uses [`Levels.offsets`](@ref) relative to that pair.

`coupling` is a complex angular Rabi-frequency matrix over `basis(s)`, with
entries only in the upper-row, lower-column block. It may come from
`Levels.coupling_matrix(s, …)`, from a relative matrix scaled with
`Levels.rabi_normalised`, or from manually specified eigenstate couplings.
A coupled-basis operator must first pass through `Levels.to_eigenbasis`.
`drives` contains [`HarmonicDrive`](@ref)s at `drive_frequency`, for example
from [`zeeman_drives`](@ref)`(s, b)`. All input matrices use the structure's
field eigenbasis; the constructor performs no rotation.

The basis must consist of exactly the two complete fine-structure manifolds
of the probe, each contiguous. Drives must be block-diagonal in these
manifolds. Complete manifolds are required for the within-manifold Floquet
dressing; see [`Levels.iscomplete`](@ref).
"""
function DrivenTransition(
    s::Eigenstructure,
    probe::Pair;
    drive_frequency,
    coupling,
    drives=HarmonicDrive[],
)
    basis = s.basis
    lower_state = validate_state_kind(s.species, probe.first)
    upper_state = validate_state_kind(s.species, probe.second)
    fs_lower = fine_structure(lower_state.level)
    fs_upper = fine_structure(upper_state.level)
    if fs_lower == fs_upper
        throw(
            ArgumentError(
                "Probed states must belong to two different fine-structure levels",
            ),
        )
    end
    manifolds = unique!([fine_structure(level) for level in basis.levels])
    if Set(manifolds) != Set([fs_lower, fs_upper])
        throw(
            ArgumentError(
                "The structure must consist of exactly the two probed manifolds " *
                "(spectator manifolds are not supported)",
            ),
        )
    end
    for fs in (fs_lower, fs_upper)
        iscomplete(s, fs) || throw(
            ArgumentError(
                "The structure must contain the complete '$fs' manifold: the " *
                "Floquet dressing runs over the whole manifold",
            ),
        )
    end
    lower_range = staterange(basis, fs_lower)
    upper_range = staterange(basis, fs_upper)
    validate_driven_blocks(length(basis), coupling, drives, lower_range, upper_range)

    # Rotating-frame diagonal at δ = 0: eigen-energies relative to the probed pair.
    lower = stateindex(basis, lower_state)
    upper = stateindex(basis, upper_state)
    energies = offsets(s)
    frame = [
        uconvert(u"rad/µs", energies[i] - energies[i in lower_range ? lower : upper]) for i in 1:length(basis)
    ]

    coupling, drives = normalise_driven_matrices(coupling, drives)

    DrivenTransition(;
        basis,
        drive_frequency,
        frame,
        drives,
        coupling,
        lower,
        upper,
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
