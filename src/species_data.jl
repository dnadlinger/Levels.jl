using Unitful.DefaultSymbols

"""
Converts a spectroscopic wavenumber (inverse wavelength) to the transition
photon energy.
"""
σ_to_energy(σ) = u"ħ" * 2π * u"c" * σ

"""
88Sr+ ion.
"""
const sr88 = NoHyperfineOneElectronSpecies(
    Dict(
        "S_1/2" => 0J,
        "D_3/2" => σ_to_energy(14555.90 / cm), # [Sansonetti2012]
        "D_5/2" => σ_to_energy(14836.24 / cm), # [Sansonetti2012]
        "P_1/2" => σ_to_energy(23715.19 / cm), # [Sansonetti2012]
        "P_3/2" => σ_to_energy(24516.65 / cm), # [Sansonetti2012]
    ),
    Dict(
        ("S_1/2", "P_3/2") => 141MHz, # …(2) [Sansonetti2012]
        ("S_1/2", "P_1/2") => 127.9MHz, # …(1.3) [Likforman2016]
        ("S_1/2", "D_5/2") => 2.559Hz, # …(10) [Sansonetti2012]
        ("S_1/2", "D_3/2") => 2.299Hz, # …(21) [Sansonetti2012]
        ("D_3/2", "P_1/2") => 7.46MHz, # …(14) [Likforman2016]
        ("D_3/2", "P_3/2") => 1.0MHz, # …(2) [Sansonetti2012]
        ("D_5/2", "P_3/2") => 8.7MHz, # …(15) [Sansonetti2012]
    )
)

export sr88
