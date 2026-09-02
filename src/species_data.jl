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
- `[Sansonetti2010]`: J. E. Sansonetti and G. Nave, "Wavelengths, Transition
  Probabilities, and Energy Levels for the Spectrum of Neutral Strontium (Sr I)",
  J. Phys. Chem. Ref. Data **39**, 033103 (2010),
  [doi:10.1063/1.3449176](https://doi.org/10.1063/1.3449176).
- `[Jiang2009]`: D. Jiang, B. Arora, M. S. Safronova, and C. W. Clark,
  "Blackbody-radiation shift in a ⁸⁸Sr⁺ ion optical frequency standard", J. Phys. B
  **42**, 154020 (2009),
  [doi:10.1088/0953-4075/42/15/154020](https://doi.org/10.1088/0953-4075/42/15/154020).
- `[Muraz2026]`: B. Muraz, M. Pepin, C. Guimard, M. Brune, B. Bakkali-Hassani,
  and S. Gleyzes, "Measuring the Sr⁺ 5s₁/₂ Landé g-factor Using Singlet-Triplet
  Oscillations in a Circular Rydberg State of Strontium",
  [arXiv:2608.26784](https://arxiv.org/abs/2608.26784) (2026).
- `[Barwood2012]`: the ⁸⁸Sr⁺ S₁/₂/D₅/₂ g-factor ratio 1.668057 reported by NPL
  at CPEM 2012 (conference presentation; no published uncertainty).
- `[Hoffman2013]`: M. R. Hoffman, T. W. Noel, C. Auchter, A. Jayakumar,
  S. R. Williams, B. B. Blinov, and E. N. Fortson, "Radio frequency spectroscopy
  measurement of the Landé g factor of the 5D₅/₂ state of Ba⁺ with a single trapped
  ion", Phys. Rev. A **88**, 025401 (2013),
  [doi:10.1103/PhysRevA.88.025401](https://doi.org/10.1103/PhysRevA.88.025401).
- `[Marx1998]`: G. Marx, G. Tommaseo, and G. Werth, "Precise g_J- and g_I-factor
  measurements of Ba⁺ isotopes", Eur. Phys. J. D **4**, 279 (1998),
  [doi:10.1007/s100530050210](https://doi.org/10.1007/s100530050210).
"""
const sr88 = NoHyperfineOneElectronSpecies(;
    # Mass of the actual ion: the neutral-atom mass 87.905612253(6) u [AME2020]
    # less one electron, plus the Sr I first ionisation energy 45932.2036(10) cm⁻¹
    # [Sansonetti2010] as its mass-equivalent binding correction — at 6.1e-9 u the
    # same size as the [AME2020] mass uncertainty.
    mass=uconvert(
        u"u",
        87.905612253u"u" - Unitful.me + σ_to_energy(45932.2036 / u"cm") / u"c"^2,
    ),
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
            # Given as the matrix element: the all-order ⟨4d₅/₂‖d‖5p₃/₂⟩ =
            # 4.187 e a₀ of [Jiang2009] (equivalent to A = 8.06 µs⁻¹) is far
            # better determined than the 8.7(15) µs⁻¹ compilation value of
            # [Sansonetti2012] (with which it is consistent, at 0.4σ), being
            # corroborated at the percent level — through its dominant share
            # of α₀(4d₅/₂) — by the measured static differential
            # polarisability of the clock transition (see
            # test-polarisability.jl).
            ("D_5/2", "P_3/2") => ReducedDipole(4.187DIPOLE_AU),
        ]
    ),
    lande_g_overrides=Dict(
        convert(NoHyperfineNumberSpec, k) => v for (k, v) in [
            "S_1/2" => 2.002290, # …(1) [Muraz2026]
            # As of 2026-09, there is no good published value for the D_5/2 g-factor
            # (or, alternatively, the S_1/2 to D_5/2 ratio). Instead, we try to get a
            # better value by naive extrapolation from ⁴⁰Ca⁺: There, the measured
            # S_1/2 / D_5/2 g-factor ratio [McMahon2026] exceeds its lowest-order
            # LS-coupling value (reduced-mass-corrected g_L, free-electron g_s, cf.
            # ls_lande_g, with the ⁴⁰Ca⁺ mass) 1.6679699 by 1.1791e-4; applying the
            # same offset to the ⁸⁸Sr⁺ LS-coupling ratio 1.6679616 gives 1.6680795,
            # i.e. g_D = 2.002290 / 1.6680795 = 1.200356 with the [Muraz2026] g_S.
            # This estimate is probably accurate to a few 1e-5. (The ratio 1.668057
            # reported by NPL at CPEM 2012 [Barwood2012] would instead give 1.2003727,
            # but no uncertainty was reported, and the number may in fact have just
            # been 1 / 0.5995. It also does not fit experimental data nearly as well.)
            "D_5/2" => 1.200356,
        ]
    ),
    polarisabilities=Dict(
        convert(NoHyperfineNumberSpec, k) => v for (k, v) in [
            # The 4d₅/₂ → 4f matrix elements and the static totals below are
            # the relativistic all-order values of [Jiang2009], Tables 1 and
            # 3; the S/D → P channel dipoles are derived from the einstein_as
            # above (cf. ImplicitPolarisability — the measured rates for the
            # S → P channels agree with the [Jiang2009] dipoles to ≈0.1%,
            # and D_5/2 → P_3/2 is [Jiang2009]-sourced there), and the
            # lumped static remainders follow by subtracting the explicit
            # channels from the totals. Those remainders cover every
            # contribution of that paper not explicit here — dominated by
            # the ionic core (5.81) plus, for D_5/2, the 5f…12f tail —
            # evaluated in the static limit. That is a good approximation
            # because they all involve transitions far above the lasers of
            # interest: the closest, 4d–5f at 71066 cm⁻¹, is enhanced by
            # only 5% at 674 nm, against 12% for the explicitly kept 4d–4f.
            "S_1/2" => ImplicitPolarisability(
                ["P_1/2", "P_3/2"];
                static_scalar=91.30POLARISABILITY_AU, # [Jiang2009]
            ),
            "D_5/2" => ImplicitPolarisability(
                [
                    "P_3/2",
                    # Explicit dipoles: these channels stay out of
                    # einstein_as, where the (dominant) 4f₅/₂ → 4d₃/₂ decay
                    # is not tabulated, so partial entries would corrupt the
                    # 4f level lifetimes.
                    "F_5/2" => 0.789DIPOLE_AU, # [Jiang2009]
                    "F_7/2" => 3.528DIPOLE_AU, # [Jiang2009]
                ];
                static_scalar=62.0POLARISABILITY_AU, # [Jiang2009]
                static_tensor=-47.7POLARISABILITY_AU, # [Jiang2009]
            ),
        ]
    ),
)

"""
⁴³Ca⁺ ion (nuclear spin ``I = 7/2``).

Hyperfine centroids are referenced to the S``_{1/2}`` centroid; all hyperfine
``A``/``B`` constants are entered as stated in the cited measurements (the signs follow
from ``μ_I < 0``). The electronic g-factors and Einstein A coefficients marked as such
are ⁴⁰Ca⁺ measurements, entered uncorrected. For the Einstein A coefficients and
g(S``_{1/2}``) the isotope dependence is far below the quoted uncertainties; for
g(D``_{5/2}``), known to 5 × 10⁻⁸ in ⁴⁰Ca⁺, the leading-order reduced-mass shift (cf.
[`Levels.ls_lande_g`](@ref)) is not, but is deliberately not applied: the many-electron
recoil corrections are unknown at that level, so will have to be measured together with
the hyperfine constants anyway.

# References

- `[AME2020]`: M. Wang, W. J. Huang, F. G. Kondev, G. Audi, and S. Naimi, "The
  AME 2020 atomic mass evaluation (II). Tables, graphs and references", Chin. Phys.
  C **45**, 030003 (2021),
  [doi:10.1088/1674-1137/abddaf](https://doi.org/10.1088/1674-1137/abddaf).
- `[Pak2022]`: C. Pak, M. J. Schlitters, and S. D. Bergeson, "Improved ionization
  potential of calcium using frequency-comb-based Rydberg spectroscopy", Phys. Rev. A
  **106**, 062818 (2022),
  [doi:10.1103/PhysRevA.106.062818](https://doi.org/10.1103/PhysRevA.106.062818).
- `[Kramida2020]`: A. Kramida, "Isotope shifts in neutral and singly-ionized
  calcium", At. Data Nucl. Data Tables **133–134**, 101322 (2020),
  [doi:10.1016/j.adt.2019.101322](https://doi.org/10.1016/j.adt.2019.101322).
- `[Arbes1994]`: F. Arbes, M. Benzing, Th. Gudjons, F. Kurth, and G. Werth,
  "Precise determination of the ground state hyperfine structure splitting of
  ⁴³Ca II", Z. Phys. D **31**, 27 (1994),
  [doi:10.1007/BF01426573](https://doi.org/10.1007/BF01426573).
- `[Nortershauser1998]`: W. Nörtershäuser et al., "Isotope shifts and hyperfine
  structure in the 3d ²D_J → 4p ²P_J transitions in calcium II", Eur. Phys. J. D
  **2**, 33 (1998), [doi:10.1007/s100530050107](https://doi.org/10.1007/s100530050107).
- `[Benhelm2007]`: J. Benhelm, G. Kirchmair, U. Rapol, T. Körber, C. F. Roos, and
  R. Blatt, "Measurement of the hyperfine structure of the S₁/₂–D₅/₂ transition in
  ⁴³Ca⁺", Phys. Rev. A **75**, 032506 (2007),
  [doi:10.1103/PhysRevA.75.032506](https://doi.org/10.1103/PhysRevA.75.032506); the
  signs of the D``_{5/2}`` constants per the erratum, Phys. Rev. A **75**, 049901
  (2007), [doi:10.1103/PhysRevA.75.049901](https://doi.org/10.1103/PhysRevA.75.049901).
- `[Tommaseo2003]`: G. Tommaseo, T. Pfeil, G. Revalde, G. Werth, P. Indelicato,
  and J. P. Desclaux, "The g_J-factor in the ground state of Ca⁺", Eur. Phys. J. D
  **25**, 113 (2003), [doi:10.1140/epjd/e2003-00096-6](https://doi.org/10.1140/epjd/e2003-00096-6).
- `[Chwalla2009]`: M. Chwalla et al., "Absolute frequency measurement of the
  ⁴⁰Ca⁺ 4s ²S₁/₂ – 3d ²D₅/₂ clock transition", Phys. Rev. Lett. **102**, 023002
  (2009), [doi:10.1103/PhysRevLett.102.023002](https://doi.org/10.1103/PhysRevLett.102.023002).
- `[Ma2024]`: Z. Ma, B. Zhang, Y. Huang, R. Hu, M. Zeng, K. Gao, and H. Guan,
  "Precision determination of the oscillating-magnetic-field-induced second-order
  Zeeman shift of a single-⁴⁰Ca⁺-ion optical clock", Phys. Rev. A **110**, 063102
  (2024), [doi:10.1103/PhysRevA.110.063102](https://doi.org/10.1103/PhysRevA.110.063102).
- `[Zhang2026]`: B. Zhang et al., "Liquid-Nitrogen-Cooled ⁴⁰Ca⁺ Ion Optical Clock
  with a Systematic Uncertainty of 4.4 × 10⁻¹⁹", Phys. Rev. Lett. **136**, 053202
  (2026), [doi:10.1103/vngc-c1xv](https://doi.org/10.1103/vngc-c1xv).
- `[McMahon2026]`: B. J. McMahon, V. S. Sandhu, J. M. Gray, C. D. Herold,
  K. R. Brown, and B. C. Sawyer, "Dual-Platform Precision Measurement of the
  3²D₅/₂ to 4²S₁/₂ g-Factor Ratio in ⁴⁰Ca⁺",
  [arXiv:2607.07929](https://arxiv.org/abs/2607.07929) (2026).
- `[Hanley2021]`: R. K. Hanley, D. T. C. Allcock, T. P. Harty, M. A. Sepiol, and
  D. M. Lucas, "Precision measurement of the ⁴³Ca⁺ nuclear magnetic moment",
  Phys. Rev. A **104**, 052804 (2021),
  [doi:10.1103/PhysRevA.104.052804](https://doi.org/10.1103/PhysRevA.104.052804).
- `[Shao2017]`: H. Shao, Y. Huang, H. Guan, C. Li, T. Shi, and K. Gao, "Precise
  determination of the quadrupole transition matrix element of ⁴⁰Ca⁺ via
  branching-fraction and lifetime measurements", Phys. Rev. A **95**, 053415
  (2017), [doi:10.1103/PhysRevA.95.053415](https://doi.org/10.1103/PhysRevA.95.053415).
- `[Hettrich2015]`: M. Hettrich et al., "Measurement of dipole matrix elements
  with a single trapped ion", Phys. Rev. Lett. **115**, 143003 (2015),
  [doi:10.1103/PhysRevLett.115.143003](https://doi.org/10.1103/PhysRevLett.115.143003).
- `[Ramm2013]`: M. Ramm, T. Pruttivarasin, M. Kokish, I. Talukdar, and
  H. Häffner, "Precision measurement method for branching fractions of excited
  P₁/₂ states applied to ⁴⁰Ca⁺", Phys. Rev. Lett. **111**, 023004 (2013),
  [doi:10.1103/PhysRevLett.111.023004](https://doi.org/10.1103/PhysRevLett.111.023004).
- `[Gerritsma2008]`: R. Gerritsma, G. Kirchmair, F. Zähringer, J. Benhelm,
  R. Blatt, and C. F. Roos, "Precision measurement of the branching fractions of
  the 4p ²P₃/₂ decay of Ca II", Eur. Phys. J. D **50**, 13 (2008),
  [doi:10.1140/epjd/e2008-00196-9](https://doi.org/10.1140/epjd/e2008-00196-9).
- `[Meir2020]`: Z. Meir, M. Sinhal, M. S. Safronova, and S. Willitsch, "Combining
  experiments and relativistic theory for establishing accurate radiative
  quantities in atoms: The lifetime of the ²P₃/₂ state in ⁴⁰Ca⁺", Phys. Rev. A
  **101**, 012509 (2020),
  [doi:10.1103/PhysRevA.101.012509](https://doi.org/10.1103/PhysRevA.101.012509).
- `[Kreuter2005]`: A. Kreuter et al., "Experimental and theoretical study of the
  3d ²D-level lifetimes of ⁴⁰Ca⁺", Phys. Rev. A **71**, 032504 (2005),
  [doi:10.1103/PhysRevA.71.032504](https://doi.org/10.1103/PhysRevA.71.032504).
- `[YuSahoo2025]`: Y. M. Yu and B. K. Sahoo, "Application of general-order
  relativistic coupled-cluster theory to estimate electric-field-response clock
  properties of Ca⁺ and Yb⁺", Phys. Rev. A **111**, 032801 (2025),
  [doi:10.1103/PhysRevA.111.032801](https://doi.org/10.1103/PhysRevA.111.032801).
- `[Tang2013]`: Y.-B. Tang, H.-X. Qiao, T.-Y. Shi, and J. Mitroy, "Dynamic
  polarizabilities for the low-lying states of Ca⁺", Phys. Rev. A **87**, 042517
  (2013),
  [doi:10.1103/PhysRevA.87.042517](https://doi.org/10.1103/PhysRevA.87.042517).
"""
const ca43 = HyperfineOneElectronSpecies(;
    # Mass of the actual ion: the neutral-atom mass 42.95876638(24) u [AME2020]
    # less one electron, plus the Ca I first ionisation energy 49305.919611(4) cm⁻¹
    # [Pak2022] as its mass-equivalent binding correction (NB: as this is well below the
    # [AME2020] mass uncertainty, the usefulness of this questionable already, so the
    # isotope shifts on the ionisation energy are especially negligible here).
    mass=uconvert(
        u"u",
        42.95876638u"u" - Unitful.me + σ_to_energy(49305.919611 / u"cm") / u"c"^2,
    ),
    nuclear_spin=7//2,
    # μ_I/μ_N = −1.315350(9)(1), the effective moment of the nucleus bound in the
    # ion, i.e. *not* corrected for diamagnetic shielding — the appropriate value
    # for the Zeeman Hamiltonian [Hanley2021]. (Older work, e.g. the [Benhelm2007]
    # fit, instead used the shielding-corrected free-nucleus value −1.317643 of
    # N. J. Stone, At. Data Nucl. Data Tables 90, 75 (2005), 0.17% away.) In the
    # H_Z = μ_B B (g_J m_J + g_I m_I) convention used here,
    # g_I = −(μ_I/μ_N)(mₑ/mₚ)/I ≈ +2.0467e-4.
    nuclear_g=(-(-1.315350) * uconvert(NoUnits, Unitful.me / Unitful.mp) / (7 // 2)),
    energies=Dict(
        convert(NoHyperfineNumberSpec, k) => v for (k, v) in [
            "S_1/2" => 0u"J",
            # Hyperfine-centroid transition frequencies from S_1/2, [Kramida2020]
            # Table 13.
            "D_3/2" => u"h" * 409_226_671.03u"MHz", # …(5), 733 nm [Kramida2020]
            "D_5/2" => u"h" * 411_046_264.4881u"MHz", # …(4), 729 nm [Kramida2020]
            "P_1/2" => u"h" * 755_223_443.81u"MHz", # …(7), 397 nm [Kramida2020]
            "P_3/2" => u"h" * 761_905_691.40u"MHz", # …(9), 393 nm [Kramida2020]
        ]
    ),
    hyperfine=Dict(
        convert(NoHyperfineNumberSpec, k) => v for (k, v) in [
            # Measured ground-state splitting Δν = 3 225 608 286.4(3) Hz
            # [Arbes1994]; A = −Δν/(I + 1/2), the sign following from μ_I < 0.
            "S_1/2" => HyperfineConstants(; a=u"h" * (-3225.60828640u"MHz") / 4),
            "P_1/2" => HyperfineConstants(; a=u"h" * -145.4u"MHz"), # …(0.1) [Nortershauser1998]
            "P_3/2" => HyperfineConstants(;
                # Note −31.0(0.2), not the −31.4 MHz sometimes transcribed
                # (e.g. in the Oxford atomic_physics package).
                a=u"h" * -31.0u"MHz", # …(0.2) [Nortershauser1998]
                b=u"h" * -6.9u"MHz", # …(1.7) [Nortershauser1998]
            ),
            "D_3/2" => HyperfineConstants(;
                a=u"h" * -47.3u"MHz", # …(0.2) [Nortershauser1998]
                b=u"h" * -3.7u"MHz", # …(1.9) [Nortershauser1998]
            ),
            # Signs per the [Benhelm2007] erratum.
            "D_5/2" => HyperfineConstants(;
                a=u"h" * -3.8931u"MHz", # …(2) [Benhelm2007]
                b=u"h" * -4.241u"MHz", # …(4) [Benhelm2007]
            ),
        ]
    ),
    lande_g_overrides=let g_s = 2.00225664 # …(9), measured in ⁴⁰Ca⁺ [Tommaseo2003]
        Dict(
            convert(NoHyperfineNumberSpec, k) => v for (k, v) in [
                "S_1/2" => g_s,
                # From the ratio g(D_5/2)/g(S_1/2) = 0.599 488 813 3(2) measured in
                # ⁴⁰Ca⁺ [McMahon2026] (Penning trap; their rf-trap value
                # 0.599 488 813(6) and the 0.599 488 79(2) of [Ma2024] concur, as does
                # the 0.599 488 818(9) of the [Zhang2026] clock evaluation — from its
                # Supplemental Material, as the raw ratio at the operating conditions,
                # whose 3e-8 excess over [Ma2024] is attributed to the trap-rf ac
                # field), i.e. 1.200 330 46(5), the uncertainty being that of g_s. This
                # supersedes the 1.2003340(3) of [Chwalla2009], whose ratio
                # 0.599 490 58(15) lies more than 10σ from all three recent
                # measurements.
                "D_5/2" => 0.5994888133 * g_s,
            ]
        )
    end,
    einstein_as=Dict(
        convert(Tuple{NoHyperfineNumberSpec,NoHyperfineNumberSpec}, k) => v for
        (k, v) in [
            # All lifetimes and branching fractions measured in ⁴⁰Ca⁺; the
            # isotope dependence is far below the quoted uncertainties.
            #
            # τ(D_5/2) = 1.1649(44) s [Shao2017]; the E2 decay to S_1/2 is the
            # only relevant channel (D_5/2 → D_3/2 M1 is ~µHz).
            ("S_1/2", "D_5/2") => 1 / 1.1649u"s",
            # τ(D_3/2) = 1176(11) ms [Kreuter2005].
            ("S_1/2", "D_3/2") => 1 / 1176u"ms",
            # τ(P_1/2) = 6.904(26) ns [Hettrich2015], split by the branching
            # fractions 0.93565(7)/0.06435(7) of [Ramm2013].
            ("S_1/2", "P_1/2") => 0.93565 / 6.904u"ns",
            ("D_3/2", "P_1/2") => 0.06435 / 6.904u"ns",
            # τ(P_3/2) = 6.639(42) ns [Meir2020] (in 6σ tension with the older
            # 6.924(19) ns of Jin & Church 1993, which [Meir2020] argues to be
            # superseded), split by the branching fractions 0.9347(3)/0.0587(2)/
            # 0.00661(4) of [Gerritsma2008].
            ("S_1/2", "P_3/2") => 0.9347 / 6.639u"ns",
            ("D_5/2", "P_3/2") => 0.0587 / 6.639u"ns",
            ("D_3/2", "P_3/2") => 0.00661 / 6.639u"ns",
        ]
    ),
    polarisabilities=Dict(
        convert(NoHyperfineNumberSpec, k) => v for (k, v) in [
            # Unlike for ⁸⁸Sr⁺, no independent sum-over-states table is
            # entered here. The explicit channel dipoles are derived (by the
            # species constructor, cf. ImplicitPolarisability) from the
            # measured Einstein A coefficients above — i.e. from
            # [Hettrich2015]/[Ramm2013] for the decays from P_1/2 and
            # [Meir2020]/[Gerritsma2008] for those from P_3/2; the derived
            # ⟨S_1/2‖d‖P_1/2⟩ = 2.8927 e a₀ equals the directly measured
            # 2.8928(43) of [Hettrich2015] — so background and near-resonant
            # channels share one source. The static totals given below are
            # the RCC values of [YuSahoo2025], Table I, whose differential
            # α₀(D_5/2) − α₀(S_1/2) = −44.02(47) agrees with the measured
            # −44.07(1) (Huang 2019); both are pinned down in
            # test-polarisability.jl.
            "S_1/2" => ImplicitPolarisability(
                ["P_1/2", "P_3/2"];
                # The implied remainder (1.58 a.u.) is less than the ionic
                # core alone (3.26, plus ≈0.19 of 5p/tail; [Tang2013] Table
                # VII): the [YuSahoo2025] total sits ≈1.9 a.u. (≈3σ) below
                # measured channels plus core, the RCC-family totals
                # generally running low against the 75.3–76.1 of the
                # DFCP/CICP/MBPT-SD/f-sum cluster ([Tang2013] Table VI). Of
                # no consequence at optical detunings, where the explicit
                # channels dominate — and largely common-mode with D_5/2
                # below, so the measured static differential stays anchored
                # regardless.
                static_scalar=74.62POLARISABILITY_AU, # …(41) [YuSahoo2025]
            ),
            "D_5/2" => ImplicitPolarisability(
                ["P_3/2"];
                # Unlike for ⁸⁸Sr⁺, the 4f/5f/higher-f channels
                # (≈54 000 cm⁻¹ up) are lumped into the implied remainder;
                # their dispersion enhancement is 7% at 729 nm and 24% at
                # 422 nm, and the DFCP decomposition puts ≈6 a.u. of static
                # f-channel content here ([Tang2013] Table X — more than the
                # whole non-core part of the remainder, cf. the S_1/2 note
                # above), bounding the lumping error at ≈0.4/1.4 a.u. at
                # those wavelengths.
                static_scalar=30.59POLARISABILITY_AU, # …(6) [YuSahoo2025]
                static_tensor=-24.50POLARISABILITY_AU, # …(12) [YuSahoo2025]
            ),
            "D_3/2" => ImplicitPolarisability(
                ["P_1/2", "P_3/2"];
                # The 3d → 4f lumping note above applies equally.
                static_scalar=33.36POLARISABILITY_AU, # …(31) [YuSahoo2025]
                static_tensor=-17.17POLARISABILITY_AU, # …(10) [YuSahoo2025]
            ),
        ]
    ),
)

export sr88, ca43
