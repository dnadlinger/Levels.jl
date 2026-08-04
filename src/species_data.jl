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

"""
⁴³Ca⁺ ion (nuclear spin ``I = 7/2``).

Hyperfine centroids are referenced to the S``_{1/2}`` centroid; all hyperfine
``A``/``B`` constants are entered as stated in the cited measurements (the signs
follow from ``μ_I < 0``). The electronic g-factors and Einstein A coefficients
marked as such are ⁴⁰Ca⁺ measurements; their isotope dependence is far below the
quoted uncertainties.

# References

- `[AME2020]`: M. Wang, W. J. Huang, F. G. Kondev, G. Audi, and S. Naimi, "The
  AME 2020 atomic mass evaluation (II). Tables, graphs and references", Chin. Phys.
  C **45**, 030003 (2021),
  [doi:10.1088/1674-1137/abddaf](https://doi.org/10.1088/1674-1137/abddaf).
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
"""
const ca43 = HyperfineOneElectronSpecies(;
    mass=42.95876638u"u", # …(24), neutral atom [AME2020]
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
    lande_g_overrides=Dict(
        convert(NoHyperfineNumberSpec, k) => v for (k, v) in [
            "S_1/2" => 2.00225664, # …(9), measured in ⁴⁰Ca⁺ [Tommaseo2003]
            "D_5/2" => 1.2003340, # …(3), measured in ⁴⁰Ca⁺ [Chwalla2009]
        ]
    ),
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
)

export sr88, ca43
