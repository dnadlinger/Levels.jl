@testitem "Gaussian beam intensity" tags = [:unit, :fast] begin
    using Unitful: mW, nm, mm, μm

    power = 1mW
    waist = 10μm
    wavelength = 400nm
    peak = 2 * power / (π * waist^2)
    rayleigh_range = π * waist^2 / wavelength

    @test gauss_intensity(power, waist, wavelength) ≈ peak
    @test gauss_intensity(power, waist, wavelength, r_offset=waist) ≈ peak * exp(-2)

    # One Rayleigh range from the focus, the beam radius has grown by √2.
    @test gauss_intensity(power, waist, wavelength, z_offset=rayleigh_range) ≈ peak / 2

    # The GaussBeamParams form defaults to the on-axis intensity at the waist,
    # wherever along the propagation direction that happens to be.
    beam = GaussBeamParams(waist, 5mm, wavelength)
    @test gauss_intensity(power, beam) ≈ peak
    @test gauss_intensity(power, beam, r_offset=waist) ≈ peak * exp(-2)
    @test gauss_intensity(power, beam, z_pos=beam.waist_pos + rayleigh_range) ≈ peak / 2

    # For a waist at the origin, z_pos and z_offset coincide.
    @test gauss_intensity(power, GaussBeamParams(waist, 0μm, wavelength), z_pos=3μm) ≈
          gauss_intensity(power, waist, wavelength, z_offset=3μm)
end
