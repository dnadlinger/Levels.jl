# Perturbative dressed-state engine: the harmonic drives dress each level
# manifold exactly (Floquet), the weak laser coupling is then treated to first
# order between the dressed Floquet modes.

"""
    dress_manifold(H0::Diagonal, H1::AbstractMatrix, Ω, harmonic_cutoff, num_time_samples;
                   min_overlap = 0.0) -> (ε, U)

Computes the Floquet modes of one manifold driven by
``H(t) = H_0 + H_1 e^{iΩt} + H_1^† e^{−iΩt}``, truncating the Floquet
Hamiltonian at `harmonic_cutoff` harmonics.

For each bare level α, the returned mode is the Floquet eigenstate with maximum
harmonic-0 weight on ``|α⟩`` (unambiguous for weak dressing); `ε[α]` is its
quasienergy (in angular frequency units, like `H0`), and `U[:, it, α]` the
periodic part ``u_α(t)`` sampled on `num_time_samples` uniform times over one drive period.

A warning is emitted whenever that maximum weight falls below `min_overlap`
(the identification of the dressed state is then ambiguous, e.g. for a
near-resonant drive).
"""
function dress_manifold(
    H0::Diagonal,
    H1::AbstractMatrix,
    Ω,
    harmonic_cutoff,
    num_time_samples;
    min_overlap=0.0,
)
    # Strip to a common inverse-time unit for the eigensolver.
    h0 = ustrip.(u"rad/µs", H0.diag)
    h1 = ustrip.(u"rad/µs", H1)
    Ω_common = ustrip(u"rad/µs", Ω)

    d = length(h0)
    dim = (2 * harmonic_cutoff + 1) * d
    F = zeros(ComplexF64, dim, dim)
    for k in (-harmonic_cutoff):harmonic_cutoff
        i0 = (k + harmonic_cutoff) * d
        for α in 1:d
            F[i0+α, i0+α] = h0[α] + k * Ω_common
        end
        if k < harmonic_cutoff
            i1 = i0 + d   # row block k + 1, column block k: coefficient of e^{+iΩt}
            F[(i1+1):(i1+d), (i0+1):(i0+d)] .= h1
            F[(i0+1):(i0+d), (i1+1):(i1+d)] .= h1'
        end
    end
    vals, vecs = eigen(Hermitian(F))

    ε = zeros(d)
    U = zeros(ComplexF64, d, num_time_samples, d)
    k0 = harmonic_cutoff * d
    θs = 2π .* (0:(num_time_samples-1)) ./ num_time_samples
    for α in 1:d
        idx = argmax(abs2.(vecs[k0+α, :]))
        overlap = abs2(vecs[k0+α, idx])
        if overlap < min_overlap
            @warn "Ambiguous dressed-state identification (e.g. near-resonant " *
                  "drive); the assignment of quasienergies to bare levels may " *
                  "be unreliable." α overlap
        end
        ε[α] = vals[idx]
        v = reshape(vecs[:, idx], d, 2 * harmonic_cutoff + 1)
        for (it, θ) in enumerate(θs), k in (-harmonic_cutoff):harmonic_cutoff
            U[:, it, α] .+= v[:, k+harmonic_cutoff+1] .* cis(k * θ)
        end
    end
    ε .* u"rad/µs", U
end

"""
Fourier components (harmonics 0 and +1) of the laser-free Hamiltonian of a
[`DrivenTransition`](@ref) restricted to the manifold covering the index range
`r`, at ``δ = 0``.
"""
manifold_blocks(dt::DrivenTransition, r::UnitRange{Int}) =
    Diagonal(dt.frame[r]), harmonic_block(dt.drives, r)

"""
A [`DrivenTransition`](@ref) with both manifolds dressed by its harmonic drives
(laser off), cached so that repeated [`sideband_rabi`](@ref) evaluations
(e.g. a modulation-index scan) reuse the diagonalisation.

The quasienergies `ε_lower`/`ε_upper` and periodic Floquet modes
`u_lower`/`u_upper` follow the conventions of [`dress_manifold`](@ref).
"""
Base.@kwdef struct DressedTransition{DT<:DrivenTransition,E<:Quantity}
    dt::DT
    ε_lower::Vector{E}
    u_lower::Array{ComplexF64,3}
    ε_upper::Vector{E}
    u_upper::Array{ComplexF64,3}
end

"""
    dress(dt::DrivenTransition, method = DressedStates())

Dresses the two manifolds of the given model by its harmonic drives, returning a
[`DressedTransition`](@ref).

`method.harmonic_cutoff` sets the Floquet truncation;
`method.num_time_samples` sets the uniform time sampling per period for the
Fourier analysis of the laser coupling in [`sideband_rabi`](@ref).
"""
function dress(dt::DrivenTransition, method::DressedStates=DressedStates())
    harmonic_cutoff = method.harmonic_cutoff
    num_time_samples = method.num_time_samples
    ε_lower, u_lower = dress_manifold(
        manifold_blocks(dt, dt.lower_range)...,
        dt.drive_frequency,
        harmonic_cutoff,
        num_time_samples,
    )
    ε_upper, u_upper = dress_manifold(
        manifold_blocks(dt, dt.upper_range)...,
        dt.drive_frequency,
        harmonic_cutoff,
        num_time_samples,
    )
    DressedTransition(; dt, ε_lower, u_lower, ε_upper, u_upper)
end

"""
    sideband_rabi(d::DressedTransition; sideband::Int, modulation = none)

Computes the first-order laser coupling between the dressed probed states at
harmonic `sideband` of the drive frequency (`+1`: probing at ν + Ω, `-1`: ν − Ω,
`0`: carrier), with the laser phase-modulated by `modulation` (a
[`HarmonicPhaseModulation`](@ref), or any callable mapping the unitful time to
the phase in radians).

Returns a NamedTuple with the effective Rabi frequency `Ω = 2|g_n|`, the
resonance detuning `δ_res` (including the ac shifts of the dressed levels), the
complex amplitude `amp = g_n` (for interference diagnostics), and the carrier
Rabi frequency `carrier`.
"""
function sideband_rabi(
    d::DressedTransition;
    sideband::Integer,
    modulation=HarmonicPhaseModulation(),
)
    dt = d.dt
    α_lower = dt.lower - first(dt.lower_range) + 1
    α_upper = dt.upper - first(dt.upper_range) + 1
    L = dt.coupling[dt.upper_range, dt.lower_range]

    num_time_samples = size(d.u_lower, 2)
    θs = 2π .* (0:(num_time_samples-1)) ./ num_time_samples
    phase_fn = phase_function(modulation, dt.drive_frequency)
    g = map(enumerate(θs)) do (it, θ)
        t = ustrip(u"s", θ * u"rad" / dt.drive_frequency)
        f = 0.5 * cis(-phase_fn(t))
        f * dot(view(d.u_upper, :, it, α_upper), L, view(d.u_lower, :, it, α_lower))
    end
    gn(n) =
        sum(g[it] * cis(-n * θs[it]) for it in 1:num_time_samples) / num_time_samples

    a = gn(sideband)
    (
        Ω=uconvert(u"rad/µs", 2 * abs(a)),
        δ_res=uconvert(
            u"rad/µs",
            d.ε_upper[α_upper] - d.ε_lower[α_lower] + sideband * dt.drive_frequency,
        ),
        amp=uconvert(u"rad/µs", a),
        carrier=uconvert(u"rad/µs", 2 * abs(gn(0))),
    )
end

"""
    sideband_rabi(dt::DrivenTransition, method = DressedStates(); sideband, modulation = none)

Convenience wrapper: dresses the transition and evaluates
[`sideband_rabi`](@ref).
"""
function sideband_rabi(
    dt::DrivenTransition,
    method::DressedStates=DressedStates();
    sideband::Integer,
    modulation=HarmonicPhaseModulation(),
)
    sideband_rabi(dress(dt, method); sideband, modulation)
end

export DressedTransition, dress, sideband_rabi
public dress_manifold
