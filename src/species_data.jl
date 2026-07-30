using Unitful

"""
Converts a spectroscopic wavenumber (inverse wavelength) to the transition
photon energy.
"""
σ_to_energy(σ) = u"ħ" * 2π * u"c" * σ

"""
⁸⁸Sr⁺ ion.

# References

- `[Sansonetti2012]`: J. E. Sansonetti, "Wavelengths, Transition Probabilities, and
  Energy Levels for the Spectra of Strontium Ions (Sr II through Sr XXXVIII)",
  J. Phys. Chem. Ref. Data **41**, 013102 (2012),
  [doi:10.1063/1.3659413](https://doi.org/10.1063/1.3659413).
- `[Likforman2016]`: J.-P. Likforman, V. Tugayé, S. Guibal, and L. Guidoni,
  "Precision measurement of the branching fractions of the 5p ²P₁/₂ state in ⁸⁸Sr⁺
  with a single ion in a microfabricated surface trap", Phys. Rev. A **93**, 052507
  (2016), [doi:10.1103/PhysRevA.93.052507](https://doi.org/10.1103/PhysRevA.93.052507).
- `[AME2020]`: M. Wang, W. J. Huang, F. G. Kondev, G. Audi, and S. Naimi, "The
  AME 2020 atomic mass evaluation (II). Tables, graphs and references", Chin. Phys.
  C **45**, 030003 (2021),
  [doi:10.1088/1674-1137/abddaf](https://doi.org/10.1088/1674-1137/abddaf).
- `[Jiang2009]`: D. Jiang, B. Arora, M. S. Safronova, and C. W. Clark,
  "Blackbody-radiation shift in a ⁸⁸Sr⁺ ion optical frequency standard", J. Phys. B
  **42**, 154020 (2009),
  [doi:10.1088/0953-4075/42/15/154020](https://doi.org/10.1088/0953-4075/42/15/154020).
"""
const sr88 = NoHyperfineOneElectronSpecies(;
    mass=87.905612253u"u", # …(6), neutral atom [AME2020]
    energies=Dict(
        convert(NoHyperfineNumberSpec, k) => v for (k, v) in [
            "S_1/2" => 0u"J",
            "D_3/2" => σ_to_energy(14555.90 / u"cm"), # [Sansonetti2012]
            "D_5/2" => σ_to_energy(14836.24 / u"cm"), # [Sansonetti2012]
            "P_1/2" => σ_to_energy(23715.19 / u"cm"), # [Sansonetti2012]
            "P_3/2" => σ_to_energy(24516.65 / u"cm"), # [Sansonetti2012]
            # 4f; only relevant as an intermediate level for the D_5/2 light shift.
            # The fine structure is inverted, F_7/2 lying below F_5/2.
            "F_5/2" => σ_to_energy(60991.34 / u"cm"), # [Sansonetti2012]
            "F_7/2" => σ_to_energy(60990.04 / u"cm"), # [Sansonetti2012]
        ]
    ),
    einstein_as=Dict(
        convert(Tuple{NoHyperfineNumberSpec,NoHyperfineNumberSpec}, k) => v for
        (k, v) in [
            ("S_1/2", "P_3/2") => 141u"µs^-1", # …(2) [Sansonetti2012]
            ("S_1/2", "P_1/2") => 127.9u"µs^-1", # …(1.3) [Likforman2016]
            ("S_1/2", "D_5/2") => 2.559u"s^-1", # …(10) [Sansonetti2012]
            ("S_1/2", "D_3/2") => 2.299u"s^-1", # …(21) [Sansonetti2012]
            ("D_3/2", "P_1/2") => 7.46u"µs^-1", # …(14) [Likforman2016]
            ("D_3/2", "P_3/2") => 1.0u"µs^-1", # …(2) [Sansonetti2012]
            ("D_5/2", "P_3/2") => 8.7u"µs^-1", # …(15) [Sansonetti2012]
        ]
    ),
    polarisabilities=Dict(
        convert(NoHyperfineNumberSpec, k) => v for (k, v) in [
            # All reduced matrix elements and static remainders below are the
            # relativistic all-order values of [Jiang2009], Tables 1 and 3.
            #
            # The remainders lump together every contribution that paper lists but
            # that is not given explicitly here, evaluated in the static limit. That
            # is a good approximation because they all involve transitions far above
            # the lasers of interest: the closest, 4d–5f at 71066 cm⁻¹, is enhanced
            # by only 5% at 674 nm, against 12% for the explicitly kept 4d–4f.
            "S_1/2" => LevelPolarisability(
                ["P_1/2" => 3.078DIPOLE_AU, "P_3/2" => 4.351DIPOLE_AU];
                # Core 5.81 + valence-core −0.26 + tail 0.02 + 6p…8p 0.016.
                static_scalar=5.586POLARIZABILITY_AU,
            ),
            "D_5/2" => LevelPolarisability(
                [
                    "P_3/2" => 4.187DIPOLE_AU,
                    "F_5/2" => 0.789DIPOLE_AU,
                    "F_7/2" => 3.528DIPOLE_AU,
                ];
                # Core 5.81 + valence-core −0.40 + tail 2.06 + 6p…8p 0.016
                # + 5f…12f 3.480.
                static_scalar=10.966POLARIZABILITY_AU,
                # Tail −0.59 + 6p…8p −0.016 + 5f…12f −0.996.
                static_tensor=-1.602POLARIZABILITY_AU,
            ),
        ]
    ),
)

export sr88
