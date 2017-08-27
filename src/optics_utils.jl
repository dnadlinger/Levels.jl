using Unitful

type GaussBeamParams
    waist_radius
    waist_pos
    wavelength
end

"""
Returns the additive unit (zero) of the same type as the given value,
to transparently support unitful and symbolic use cases.
"""
zero_like(x) = x - x

"""
Returns the intensity at the given point of a Gaussian beam with the
specified parameters and total power.

Parameters
----------
power
    The optical power integrated over the whole beam profile.
waist
    The :math:`1 / e^2`-radius of the beam at its waist.
wavelength
    The wavelength of the radiation.
r_offset
    If given, the distance of the point from the focus in the radial direction.
z_offset
    If given, the distance of the point from the focus along the direction of
    propagation.
"""
function gauss_intensity(power, waist, wavelength; r_offset=zero_like(waist),
        z_offset=zero_like(waist))
    gauss_intensity(power, GaussBeamParams(waist, zero_like(waist), wavelength),
        r_offset=r_offset, z_pos=z_offset)
end

function gauss_intensity(power, beam::GaussBeamParams; r_offset=zero_like(waist),
        z_pos=zero_like(waist))
    raleigh_range = π * beam.waist_radius^2 / beam.wavelength
    radius = beam.waist_radius * sqrt(1 + ((z_pos - beam.waist_pos) / raleigh_range)^2)
    2 * power / (π * radius^2) * exp(-2 *(r_offset / radius)^2)
end

export GaussBeamParams, gauss_intensity
