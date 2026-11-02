# Specification types: the laser beams, the laser scheme (beams and rotating-frame
# choice) and the motional modes — parameters only, no atomic data; the atom is
# the Eigenstructure the scheme is applied to.

"""
One running-wave laser beam.

The frequency is a [`RelativeFrequency`](@ref) — a detuning from the zero-field
interval of a named pair of levels (``F`` levels for a hyperfine species), so
the tabulated line centres never enter the model. `ε` and `direction` are the
(complex) polarisation and the unit propagation direction in the frame with the
quantisation axis along ẑ (cf. [`Levels.beam_vectors`](@ref)); only their directions
matter, the field amplitude being fixed by `intensity`. `linewidth` is the
Lorentzian full width at half maximum of the laser spectrum (angular units),
modelled as phase diffusion (cf. [`lindblad_model`](@ref)).

Constructed via `LaserBeam(frequency, intensity, ε, direction; linewidth)`, or with the
level pair and detuning spelled out as `LaserBeam(lower => upper, detuning,
intensity, ε, direction; linewidth)`.
"""
struct LaserBeam{F<:RelativeFrequency,I<:Quantity,W<:Quantity}
    "Laser frequency relative to a zero-field line centre."
    frequency::F

    "Beam intensity at the ion (W/m²)."
    intensity::I

    "Polarisation vector (unit length, complex)."
    ε::SVector{3,ComplexF64}

    "Propagation direction (unit length)."
    direction::SVector{3,Float64}

    "Laser linewidth (FWHM, angular units)."
    linewidth::W
end

function LaserBeam(
    frequency::RelativeFrequency,
    intensity,
    ε,
    direction;
    linewidth=0.0u"rad/µs",
)
    if dimension(intensity) != dimension(1.0u"W/m^2") ||
       !isreal(intensity) ||
       !isfinite(intensity) ||
       intensity < zero(intensity)
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
                "(e.g. angular(100u\"kHz\")), got $linewidth",
            ),
        )
    end
    if length(ε) != 3 || length(direction) != 3
        throw(ArgumentError("Polarisation and beam direction must be 3-vectors"))
    end
    ε_scale = sqrt(sum(abs2, ε))
    n_scale = sqrt(sum(abs2, direction))
    if iszero(ε_scale) || iszero(n_scale)
        throw(ArgumentError("Polarisation and beam direction must be non-zero"))
    end
    if !all(isreal, direction)
        throw(ArgumentError("The beam direction must be a real vector"))
    end
    ε_unit = SVector{3,ComplexF64}(ε ./ ε_scale)
    n_unit = SVector{3,Float64}(real.(direction) ./ n_scale)
    if abs(sum(ε_unit .* n_unit)) > 1e-6
        throw(
            ArgumentError(
                "Polarisation must be transverse to the beam direction (ε ⋅ direction = 0)",
            ),
        )
    end
    LaserBeam(frequency, intensity, ε_unit, n_unit, linewidth)
end

LaserBeam(pair::Pair, detuning, intensity, ε, direction; linewidth=0.0u"rad/µs") =
    LaserBeam(RelativeFrequency(pair, detuning), intensity, ε, direction; linewidth)

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
A collection of laser beams and a rotating-frame choice, applied to an
[`Eigenstructure`](@ref Levels.Eigenstructure) by [`lindblad_model`](@ref).

Constructed with `LaserScheme(beams; frame_reference, frame_kind)`.
`frame_kind` selects one frame frequency per fine-structure level (`:levels`,
default) or per state (`:states`); see [`RotatingFrame`](@ref).
`frame_reference` sets the root level's frame frequency to its zero-field
centroid. It defaults to the upper level of the first beam, or the first
basis level if there are no beams. Changing this reference leaves observables
unchanged. Beams and reference levels are validated when the frame is built.
"""
struct LaserScheme
    beams::Vector{LaserBeam}
    frame_reference::Union{Nothing,NoHyperfineNumberSpec}
    frame_kind::Symbol
end

function LaserScheme(beams; frame_reference=nothing, frame_kind::Symbol=:levels)
    beams = LaserBeam[beams...]
    frame_kind in (:levels, :states) ||
        throw(ArgumentError("frame_kind must be :levels or :states, got $frame_kind"))
    reference = if isnothing(frame_reference)
        isempty(beams) ? nothing : beam_levels(first(beams))[2]
    else
        fine_structure(frame_reference)
    end
    LaserScheme(beams, reference, frame_kind)
end

"""
One motional mode of the trapped ion: its (angular) frequency, unit direction
in the frame with the quantisation axis along ẑ, and `participation`.

`participation` ``b`` is the magnitude of the modelled ion's component in the
mass-weighted normal-mode vector of a crystal, so that its position along
`direction` is ``b\\sqrt{ħ/2Mω}(a + a^†)`` with ``M`` the ion's own mass. It
scales every Lamb–Dicke factor of the mode, beams and recoil alike; 1 is a
single ion. A sign of the mode vector at this ion belongs in `direction`.

Constructed via `MotionalMode(frequency, direction; participation = 1)`; the
direction is normalised.
"""
struct MotionalMode{F<:Quantity}
    frequency::F
    direction::SVector{3,Float64}
    participation::Float64
end

function MotionalMode(frequency, direction; participation::Real=1)
    (isfinite(participation) && participation > 0) || throw(
        ArgumentError(
            "The mode participation must be positive and finite, got $participation",
        ),
    )
    if !(frequency isa Unitful.Frequency) || frequency <= zero(frequency)
        throw(
            ArgumentError(
                "The mode frequency must be a positive angular frequency " *
                "(e.g. angular(1.2u\"MHz\")), got $frequency",
            ),
        )
    end
    if length(direction) != 3 || !all(isreal, direction) || all(iszero, direction)
        throw(ArgumentError("The mode direction must be a non-zero real 3-vector"))
    end
    d = SVector{3,Float64}(real.(direction))
    MotionalMode(frequency, d ./ sqrt(sum(abs2, d)), Float64(participation))
end

export LaserBeam, LaserScheme, MotionalMode
