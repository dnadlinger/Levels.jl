# Perturbative dressed-state engine: the harmonic drives dress each level
# manifold exactly (Floquet), the weak laser coupling is then treated to first
# order between the dressed Floquet modes.

"""
    dress_manifold(H0::Diagonal, H1::AbstractMatrix, Ω, nharm, ngrid) -> (ε, U)

Computes the Floquet modes of one manifold driven by
``H(t) = H_0 + H_1 e^{iΩt} + H_1^† e^{−iΩt}``, truncating the Floquet
Hamiltonian at `nharm` harmonics.

For each bare level α, the returned mode is the Floquet eigenstate with maximum
harmonic-0 weight on ``|α⟩`` (unambiguous for weak dressing); `ε[α]` is its
quasienergy (in angular frequency units, like `H0`), and `U[:, it, α]` the
periodic part ``u_α(t)`` sampled on `ngrid` uniform times over one drive period.
"""
function dress_manifold(H0::Diagonal, H1::AbstractMatrix, Ω, nharm, ngrid)
    # Strip to a common inverse-time unit for the eigensolver.
    h0 = ustrip.(u"µs^-1", H0.diag)
    h1 = ustrip.(u"µs^-1", H1)
    Ω_common = ustrip(u"µs^-1", Ω)

    d = length(h0)
    dim = (2 * nharm + 1) * d
    F = zeros(ComplexF64, dim, dim)
    for k in (-nharm):nharm
        i0 = (k + nharm) * d
        for α in 1:d
            F[i0+α, i0+α] = h0[α] + k * Ω_common
        end
        if k < nharm
            i1 = i0 + d   # row block k + 1, column block k: coefficient of e^{+iΩt}
            F[(i1+1):(i1+d), (i0+1):(i0+d)] .= h1
            F[(i0+1):(i0+d), (i1+1):(i1+d)] .= h1'
        end
    end
    vals, vecs = eigen(Hermitian(F))

    ε = zeros(d)
    U = zeros(ComplexF64, d, ngrid, d)
    k0 = nharm * d
    θs = 2π .* (0:(ngrid-1)) ./ ngrid
    for α in 1:d
        idx = argmax(abs2.(vecs[k0+α, :]))
        ε[α] = vals[idx]
        v = reshape(vecs[:, idx], d, 2 * nharm + 1)
        for (it, θ) in enumerate(θs), k in (-nharm):nharm
            U[:, it, α] .+= v[:, k+nharm+1] .* cis(k * θ)
        end
    end
    ε .* u"µs^-1", U
end

"""
Fourier components (harmonics 0 and +1) of the laser-free Hamiltonian of a
[`DrivenTransition`](@ref) restricted to the manifold covering the index range
`r`, at ``δ = 0``.
"""
function manifold_blocks(dt::DrivenTransition, r::UnitRange{Int})
    H0 = Diagonal(dt.frame[r])
    H1 = zeros(ComplexF64, length(r), length(r)) .* u"µs^-1"
    for drive in dt.drives
        H1 .+= (0.5 * cis(drive.phase)) .* drive.amplitude[r, r]
    end
    H0, H1
end

"""
A [`DrivenTransition`](@ref) with both manifolds dressed by its harmonic drives
(laser off), cached so that repeated [`sideband_amplitude`](@ref) evaluations
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
    dress(dt::DrivenTransition; nharm = 8, ngrid = 256)

Dresses the two manifolds of the given model by its harmonic drives, returning a
[`DressedTransition`](@ref).

`nharm` is the Floquet harmonic truncation of the dressing; `ngrid` the number of
uniform time samples per drive period used for the Fourier analysis of the laser
coupling in [`sideband_amplitude`](@ref).
"""
function dress(dt::DrivenTransition; nharm::Int=8, ngrid::Int=256)
    ε_lower, u_lower = dress_manifold(
        manifold_blocks(dt, dt.lower_range)...,
        dt.drive_frequency,
        nharm,
        ngrid,
    )
    ε_upper, u_upper = dress_manifold(
        manifold_blocks(dt, dt.upper_range)...,
        dt.drive_frequency,
        nharm,
        ngrid,
    )
    DressedTransition(; dt, ε_lower, u_lower, ε_upper, u_upper)
end

"""
    sideband_amplitude(d::DressedTransition; sideband::Int, modulation = none)

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
function sideband_amplitude(
    d::DressedTransition;
    sideband::Int,
    modulation=HarmonicPhaseModulation(),
)
    dt = d.dt
    α_lower = dt.lower - first(dt.lower_range) + 1
    α_upper = dt.upper - first(dt.upper_range) + 1
    L = dt.coupling[dt.upper_range, dt.lower_range]

    ngrid = size(d.u_lower, 2)
    θs = 2π .* (0:(ngrid-1)) ./ ngrid
    phase_fn = phase_function(modulation, dt.drive_frequency)
    g = map(enumerate(θs)) do (it, θ)
        t = ustrip(u"s", θ / dt.drive_frequency)
        f = 0.5 * cis(-phase_fn(t))
        f * dot(view(d.u_upper, :, it, α_upper), L, view(d.u_lower, :, it, α_lower))
    end
    gn(n) = sum(g[it] * cis(-n * θs[it]) for it in 1:ngrid) / ngrid

    a = gn(sideband)
    (
        Ω=uconvert(u"ms^-1", 2 * abs(a)),
        δ_res=uconvert(
            u"µs^-1",
            d.ε_upper[α_upper] - d.ε_lower[α_lower] + sideband * dt.drive_frequency,
        ),
        amp=uconvert(u"ms^-1", a),
        carrier=uconvert(u"ms^-1", 2 * abs(gn(0))),
    )
end

"""
    sideband_rabi(dt::DrivenTransition; sideband, modulation = none, nharm = 8, ngrid = 256)

Convenience wrapper: dresses the transition and evaluates
[`sideband_amplitude`](@ref).
"""
function sideband_rabi(
    dt::DrivenTransition;
    sideband::Int,
    modulation=HarmonicPhaseModulation(),
    nharm::Int=8,
    ngrid::Int=256,
)
    sideband_amplitude(dress(dt; nharm, ngrid); sideband, modulation)
end

export DressedTransition, dress, sideband_amplitude, sideband_rabi
public dress_manifold
