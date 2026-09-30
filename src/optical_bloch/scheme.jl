# Specification types: the laser beams, the internal-state problem (basis, beams,
# static field) and the motional modes — parameters only, no atomic data.

"""
One running-wave laser beam.

The frequency is a [`RelativeFrequency`](@ref) — a detuning from the zero-field
interval of a named pair of levels (``F`` levels for a hyperfine species), so
the tabulated line centres never enter the model. `ε` and `n` are the
(complex) polarisation and the unit propagation direction in the frame with the
quantisation axis along ẑ (cf. [`Levels.beam_vectors`](@ref)); only their directions
matter, the field amplitude being fixed by `intensity`. `linewidth` is the
Lorentzian full width at half maximum of the laser spectrum (angular units),
modelled as phase diffusion (cf. [`lindblad_model`](@ref)).

Constructed via `LaserBeam(frequency, intensity, ε, n; linewidth)`, or with the
level pair and detuning spelled out as `LaserBeam(lower => upper, detuning,
intensity, ε, n; linewidth)`.
"""
struct LaserBeam{F<:RelativeFrequency,I<:Quantity,W<:Quantity}
    "Laser frequency relative to a zero-field line centre."
    frequency::F

    "Beam intensity at the ion (W/m²)."
    intensity::I

    "Polarisation vector (unit length, complex)."
    ε::SVector{3,ComplexF64}

    "Propagation direction (unit length)."
    n::SVector{3,Float64}

    "Laser linewidth (FWHM, angular units)."
    linewidth::W
end

function LaserBeam(frequency::RelativeFrequency, intensity, ε, n; linewidth=0.0u"µs^-1")
    if dimension(intensity) != dimension(1.0u"W/m^2")
        throw(
            ArgumentError(
                "The beam intensity must be a power per area, got $intensity",
            ),
        )
    end
    if !(linewidth isa Unitful.Frequency) || linewidth < zero(linewidth)
        throw(
            ArgumentError(
                "The laser linewidth must be a non-negative angular frequency " *
                "(e.g. 2π * 100u\"kHz\"), got $linewidth",
            ),
        )
    end
    if length(ε) != 3 || length(n) != 3
        throw(ArgumentError("Polarisation and beam direction must be 3-vectors"))
    end
    ε_scale = sqrt(sum(abs2, ε))
    n_scale = sqrt(sum(abs2, n))
    if iszero(ε_scale) || iszero(n_scale)
        throw(ArgumentError("Polarisation and beam direction must be non-zero"))
    end
    if !all(isreal, n)
        throw(ArgumentError("The beam direction must be a real vector"))
    end
    ε_unit = SVector{3,ComplexF64}(ε ./ ε_scale)
    n_unit = SVector{3,Float64}(real.(n) ./ n_scale)
    if abs(sum(ε_unit .* n_unit)) > 1e-6
        throw(
            ArgumentError(
                "Polarisation must be transverse to the beam direction (ε ⋅ n = 0)",
            ),
        )
    end
    LaserBeam(frequency, intensity, ε_unit, n_unit, linewidth)
end

LaserBeam(pair::Pair, detuning, intensity, ε, n; linewidth=0.0u"µs^-1") =
    LaserBeam(RelativeFrequency(pair, detuning), intensity, ε, n; linewidth)

"""
The lower and upper fine-structure levels a beam connects.
"""
beam_levels(beam::LaserBeam) =
    (fine_structure(beam.frequency.lower), fine_structure(beam.frequency.upper))

"""
Returns the fine-structure levels of a basis, in order of first appearance.
"""
fine_structure_levels(basis::StateBasis) = unique!(fine_structure.(basis.levels))

"""
The internal-state problem: the Zeeman (or hyperfine) states kept, the laser
beams driving them, and the static magnetic field along the quantisation axis
ẑ.

`frame_reference` names the level whose rotating-frame frequency is its own
zero-field centroid; the frame frequencies of the levels connected to it by
beams follow from the beam frequencies (cf. [`RotatingFrame`](@ref)). It
defaults to the upper level of the first beam. The choice does not affect any
observable — it only fixes which level's states sit at their bare Zeeman
(hyperfine) energies on the diagonal of the rotating-frame Hamiltonian.

For a hyperfine species the field must be non-zero, as the basis states denote
the adiabatically-labelled field eigenstates (cf. [`hyperfine_manifold`](@ref)).

Constructed via `LaserScheme(basis, beams; static_field, frame_reference)`; the
atomic data enter only when a [`LindbladModel`](@ref) is derived from it.
"""
struct LaserScheme{B<:StateBasis,F<:Quantity}
    basis::B
    beams::Vector{LaserBeam}
    static_field::F
    frame_reference::NoHyperfineNumberSpec
end

function LaserScheme(basis::StateBasis, beams; static_field, frame_reference=nothing)
    beams = LaserBeam[beams...]
    isempty(basis) && throw(ArgumentError("The basis must contain at least one state"))
    if dimension(static_field) != dimension(1.0u"T")
        throw(
            ArgumentError("The static field must be a flux density, got $static_field"),
        )
    end
    levels = fine_structure_levels(basis)
    level_kind = eltype(basis.levels)
    for (i, beam) in enumerate(beams)
        lo, hi = beam_levels(beam)
        for (which, level) in (("lower", lo), ("upper", hi))
            if !(level in levels)
                throw(
                    ArgumentError(
                        "Beam $i: $which level '$level' is not part of the basis",
                    ),
                )
            end
        end
        if !(beam.frequency.lower isa level_kind)
            throw(
                ArgumentError(
                    "Beam $i: the reference levels must be of the basis' kind " *
                    "('$(beam.frequency.lower)' for a basis of $level_kind levels)",
                ),
            )
        end
    end
    reference = if isnothing(frame_reference)
        isempty(beams) ? first(levels) : beam_levels(first(beams))[2]
    else
        fine_structure(frame_reference)
    end
    if !(reference in levels)
        throw(
            ArgumentError(
                "Frame reference level '$reference' is not part of the basis",
            ),
        )
    end
    LaserScheme(basis, beams, static_field, reference)
end

"""
One motional mode of the trapped ion: its (angular) frequency and unit
direction in the frame with the quantisation axis along ẑ.

Constructed via `MotionalMode(frequency, direction)`; the direction is
normalised.
"""
struct MotionalMode{F<:Quantity}
    frequency::F
    direction::SVector{3,Float64}
end

function MotionalMode(frequency, direction)
    if !(frequency isa Unitful.Frequency) || frequency <= zero(frequency)
        throw(
            ArgumentError(
                "The mode frequency must be a positive angular frequency " *
                "(e.g. 2π * 1.2u\"MHz\"), got $frequency",
            ),
        )
    end
    if length(direction) != 3 || !all(isreal, direction) || all(iszero, direction)
        throw(ArgumentError("The mode direction must be a non-zero real 3-vector"))
    end
    d = SVector{3,Float64}(real.(direction))
    MotionalMode(frequency, d ./ sqrt(sum(abs2, d)))
end

export LaserBeam, LaserScheme, MotionalMode
