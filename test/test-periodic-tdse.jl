@testitem "Floquet vs exact monodromy engine" tags=[:integration, :slow] setup=[
    PeriodicSetup,
] begin
    θ = 0.7
    drives = [
        zeeman_drive(sr88, BASIS, [0, 0, 1] .* 2.5u"µT"; phase=π / 2),
        zeeman_drive(sr88, BASIS, [cos(θ), sin(θ), 0] .* 25.0u"µT"; phase=π / 2),
    ]
    β_Q = 0.01
    for (tr, sb, β_I) in (
        (TRANSITIONS[2], +1, 0.02),
        (TRANSITIONS[7], -1, -0.005),
        (TRANSITIONS[4], -1, 0.01),
    )
        dt = make_transition(tr; drives)
        modulation = HarmonicPhaseModulation(β_I, β_Q)
        fl = sideband_rabi(dt; sideband=sb, modulation)
        ex = exact_sideband(dt; sideband=sb, modulation)
        @test fl.Ω ≈ ex.Ω rtol = 1e-3
        # δ_res agrees up to the laser-induced ac Stark shift (∝ Ω0², ~200 Hz here).
        @test abs(fl.δ_res - ex.δ_res) < 2π * 500.0u"Hz"
    end
end

@testitem "Stroboscopic Rabi flop" tags=[:integration, :slow] setup=[PeriodicSetup] begin
    θ = 0.7
    drives = [
        zeeman_drive(sr88, BASIS, [0, 0, 1] .* 2.5u"µT"; phase=π / 2),
        zeeman_drive(sr88, BASIS, [cos(θ), sin(θ), 0] .* 25.0u"µT"; phase=π / 2),
    ]
    tr = TRANSITIONS[7]
    dt = make_transition(tr; drives)
    modulation = HarmonicPhaseModulation(0.01, 0.01)

    ex = exact_sideband(dt; sideband=+1, modulation)
    period = 2π / Ω_RF
    n_π = round(Int, uconvert(NoUnits, π / ex.Ω / period))   # periods per π-pulse
    populations =
        stroboscopic_populations(dt, round(Int, 1.2 * n_π); δ=ex.δ_res, modulation)
    target = populations[stateindex(BASIS, tr.second), :]
    @test maximum(target) > 0.99                  # full transfer on resonance
    @test abs(argmax(target) - 1 - n_π) < 0.05 * n_π   # flop time matches π/Ω
end
