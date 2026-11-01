@testitem "Quasienergy shifts of a rotating two-state coupling" tags=[:unit, :fast] begin
    using LinearAlgebra
    using Levels.PeriodicDynamics
    using Unitful

    # One rotating Fourier component gives an isolated two-state Floquet
    # block, whose exact eigenvalues also fix the sideband sign convention.
    Δ = angular(5.0u"MHz")
    ω = angular(2.0u"MHz")
    g = angular(0.01u"MHz")
    H0 = Diagonal([zero(Δ), Δ])
    H1 = [0.0 im; 0.0 0.0] .* g
    detuning = Δ - ω
    pt = -g^2 / detuning
    exact = (detuning - sqrt(detuning^2 + 4g^2)) / 2

    @test quasienergy_shift(H0, H1, ω, 1) ≈ pt
    @test quasienergy_shift(H0, H1, ω, 2) ≈ -pt
    @test quasienergy_shift(H0, H1, ω, 1, Floquet(; harmonic_cutoff=2)) ≈ exact
    @test quasienergy_shift(H0, H1, ω, 2, Floquet(; harmonic_cutoff=2)) ≈ -exact
    terms = quasienergy_shift_contributions(H0, H1, ω, 1)
    @test length(terms) == 1
    @test only(terms).k == 2
    @test only(terms).sideband == 1
    @test only(terms).coupling ≈ g
    @test only(terms).detuning ≈ -detuning
    @test only(terms).shift ≈ pt
    @test unit(only(terms).shift) == u"rad/µs"
end

@testitem "Magnetic-phasor quasienergy convenience API" tags=[:unit, :fast] begin
    using Levels.PeriodicDynamics
    using Unitful

    m = Levels.hyperfine_manifold(ca43, "S_1/2"; B=0.5u"mT")
    state = StateSpec("S_1/2 F=4", 0)
    ω = angular(50u"MHz")
    for b in ([1.0, 0.0, 0.0]u"µT", [1.0, -1.0im, 0.5]u"µT")
        drives = zeeman_drives(ca43, m.basis, b)
        for method in (PerturbationTheory(), Floquet(; harmonic_cutoff=3))
            @test ac_zeeman_shift(m, state, ω, b, method) ==
                  quasienergy_shift(m, state, ω, drives, method)
        end
        @test !applicable(quasienergy_shift, m, state, ω, b)
        @test !applicable(ac_zeeman_shift, m, state, ω, drives)
    end
    @test PerturbationTheory() isa QuasienergyShiftMethod
    @test Floquet() isa QuasienergyShiftMethod
end
