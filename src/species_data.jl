using Unitful.DefaultSymbols

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
"""
const sr88 = NoHyperfineOneElectronSpecies(
    87.905612253u"u", # …(6), neutral atom [AME2020]
    Dict(
        convert(NoHyperfineNumberSpec, k) => v for (k, v) in [
            "S_1/2" => 0J,
            "D_3/2" => σ_to_energy(14555.90 / cm), # [Sansonetti2012]
            "D_5/2" => σ_to_energy(14836.24 / cm), # [Sansonetti2012]
            "P_1/2" => σ_to_energy(23715.19 / cm), # [Sansonetti2012]
            "P_3/2" => σ_to_energy(24516.65 / cm), # [Sansonetti2012]
        ]
    ),
    Dict(
        convert(Tuple{NoHyperfineNumberSpec,NoHyperfineNumberSpec}, k) => v for
        (k, v) in [
            ("S_1/2", "P_3/2") => 141MHz, # …(2) [Sansonetti2012]
            ("S_1/2", "P_1/2") => 127.9MHz, # …(1.3) [Likforman2016]
            ("S_1/2", "D_5/2") => 2.559Hz, # …(10) [Sansonetti2012]
            ("S_1/2", "D_3/2") => 2.299Hz, # …(21) [Sansonetti2012]
            ("D_3/2", "P_1/2") => 7.46MHz, # …(14) [Likforman2016]
            ("D_3/2", "P_3/2") => 1.0MHz, # …(2) [Sansonetti2012]
            ("D_5/2", "P_3/2") => 8.7MHz, # …(15) [Sansonetti2012]
        ]
    ),
)

export sr88
