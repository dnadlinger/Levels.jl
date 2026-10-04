# The motional layer: Lamb–Dicke factors, the first-order sideband Hamiltonian
# and the recoil jump operators of one motional mode, derived from a
# LindbladModel of the internal states.

"""
Coupling of one [`MotionalMode`](@ref) to the internal dynamics of a
[`LindbladModel`](@ref), to first order in the Lamb–Dicke parameters.

With ``x = x_0 (a + a^†)``, ``x_0 = \\sqrt{ħ / 2 M ω_m}``, along the mode
direction, the full Hamiltonian is ``H_\\mathrm{int} ⊗ 1 + ω_m a^† a +
H_\\mathrm{sb} ⊗ (a + a^†)`` and each spontaneous-emission jump operator ``L``
acquires the companion ``\\sqrt{α} η_0 L ⊗ (a + a^†)`` for the recoil of the
emitted photon.

Two distinct Lamb–Dicke quantities appear, and the field names spell out which:

- `projected_lamb_dicke[l]` ``= (\\vec{k}_l ⋅ \\vec{e}_m) x_0`` is the **signed,
  projected** Lamb–Dicke factor of beam ``l`` — the coefficient of ``(a + a^†)``
  in ``e^{i \\vec{k}_l ⋅ \\vec{r}}`` — which is what enters the sideband
  Hamiltonian ``H_\\mathrm{sb} = \\sum_l η_l \\, i (C_l - C_l^†)/2``. For a Raman
  process through two beams, the difference of their values is the two-photon
  factor.
- `recoil_lamb_dicke[k]` ``= |\\vec{k}| x_0`` is the **unprojected** Lamb–Dicke
  parameter of the decay transition behind jump operator ``k``; the direction of
  the emitted photon is averaged over its radiation pattern into the second
  moment ``α_k = ⟨(\\hat{k} ⋅ \\vec{e}_m)^2⟩`` (`recoil_moments`), which depends
  on the emitted spherical component ``q`` and the angle between the mode and
  the quantisation axis: ``1/5`` (``q = 0``) and ``2/5`` (``|q| = 1``) for a mode
  along ẑ, ``1/3`` on average over ``q`` (cf. [`recoil_moment`](@ref)).

`recoil_operators[k]` is ``\\sqrt{α_k} η_{0,k} L_k`` for the ``k``-th
spontaneous-emission operator of the model (`recoil_indices[k]` its index in
`model.jump_operators`; laser-dephasing operators carry no recoil). `emission`
records how the moments were obtained — `:exact` from the radiation pattern
of each emitted component, `:isotropic` for uniform emission, `:constant` for
a user-supplied number — which is what decides whether the full emission
direction distribution is available to [`motional_model`](@ref) beyond first
order in the Lamb–Dicke parameters. `sideband_harmonics` mirrors
`model.harmonics` for beams with a beat note.

Constructed via [`motional_coupling`](@ref).
"""
struct MotionalCoupling{M<:MotionalMode,H<:Quantity,W<:Quantity,J<:Quantity}
    mode::M
    projected_lamb_dicke::Vector{Float64}
    recoil_lamb_dicke::Vector{Float64}
    recoil_moments::Vector{Float64}
    emission::Symbol
    sideband_hamiltonian::Matrix{H}
    sideband_harmonics::Vector{Tuple{W,Matrix{H}}}
    recoil_operators::Vector{Matrix{J}}
    recoil_indices::Vector{Int}
end

"""
    lamb_dicke(species, mode::MotionalMode, wavevector; projected = true)
    lamb_dicke(species, mode::MotionalMode, beam::LaserBeam; projected = true)

Returns the Lamb–Dicke factor of a plane wave with the given wavevector (a
Cartesian 3-vector of inverse lengths, or the beam's ``k = ω n / c``) and the
given mode of an ion of the species' mass: the signed **projected** value
``(\\vec{k} ⋅ \\vec{e}_m) \\sqrt{ħ / 2 M ω_m}`` by default, or with
`projected = false` the **unprojected** ``|\\vec{k}| \\sqrt{ħ / 2 M ω_m}``, the
conventional Lamb–Dicke parameter of the transition (cf.
[`MotionalCoupling`](@ref) for which is which).
"""
function lamb_dicke(species, mode::MotionalMode, wavevector; projected::Bool=true)
    x0 = sqrt(u"ħ" / (2 * species.mass * mode.frequency))
    k = projected ? sum(wavevector .* mode.direction) : sqrt(sum(abs2, wavevector))
    uconvert(NoUnits, k * x0)
end

function lamb_dicke(species, mode::MotionalMode, beam::LaserBeam; projected::Bool=true)
    k = photon_energy(species, beam.frequency) / (u"ħ" * u"c")
    lamb_dicke(species, mode, k .* beam.n; projected)
end

# Angular radiation patterns of the electric multipole components, as
# polynomial coefficients in c = cos θ (θ from the quantisation axis), up to a
# constant: E1 q = 0: sin²θ, |q| = 1: (1 + cos²θ)/2; E2 q = 0: sin²θ cos²θ,
# |q| = 1: 1 − 3 cos²θ + 4 cos⁴θ, |q| = 2: 1 − cos⁴θ (Jackson, Classical
# Electrodynamics, Sec. 9.7). Rank 0 stands for isotropic emission.
function radiation_pattern(rank, q)
    if rank == 0
        (1.0,)
    elseif rank == 1
        q == 0 ? (1.0, 0.0, -1.0) : (1.0, 0.0, 1.0)
    elseif rank == 2
        if q == 0
            (0.0, 0.0, 1.0, 0.0, -1.0)
        elseif abs(q) == 1
            (1.0, 0.0, -3.0, 0.0, 4.0)
        else
            (1.0, 0.0, 0.0, 0.0, -1.0)
        end
    else
        throw(
            ArgumentError(
                "Radiation patterns are tabulated for ranks 1 and 2 (and 0 for " *
                "isotropic emission) only",
            ),
        )
    end
end

# ∫_{-1}^{1} c^n p(c) dc for a polynomial p given by its coefficients.
function moment_integral(coefficients, n)
    total = 0.0
    for (i, a) in enumerate(coefficients)
        power = n + i - 1
        iseven(power) && (total += a * 2 / (power + 1))
    end
    total
end

"""
    recoil_moment(rank, q, direction)

Returns the second moment ``α = ⟨(\\hat{k} ⋅ \\vec{e}_m)^2⟩`` of the direction of a
photon emitted in the spherical component `q` of an electric multipole
transition of the given `rank`, projected onto the unit mode `direction`
(quantisation axis along ẑ): the fraction of the squared recoil that a
spontaneous emission imparts, on average, to that mode.

By the azimuthal symmetry of the radiation patterns,
``α = ⟨\\hat{k}_z^2⟩ \\cos^2 θ_m + \\tfrac{1}{2}(1 - ⟨\\hat{k}_z^2⟩) \\sin^2 θ_m``
with ``θ_m`` the angle between mode and quantisation axis; for electric-dipole
emission ``⟨\\hat{k}_z^2⟩ = 1/5`` (``q = 0``) and ``2/5`` (``|q| = 1``), so a mode
along ẑ sees ``1/5`` and ``2/5``, one transverse to it ``2/5`` and ``3/10``,
and the average over the three components is the isotropic ``1/3`` in every
direction (cf. Javanainen & Stenholm, Appl. Phys. **21**, 35 (1980)). `rank =
0` denotes isotropic emission, with ``α = 1/3`` for every direction. The full
distribution of ``\\hat{k} ⋅ \\vec{e}_m`` behind this moment is available as a
quadrature rule from [`emission_rule`](@ref).
"""
function recoil_moment(rank, q, direction)
    abs(q) <= rank || throw(ArgumentError("|q| must not exceed the rank"))
    p = radiation_pattern(rank, q)
    kz2 = moment_integral(p, 2) / moment_integral(p, 0)
    cos2 = abs2(direction[3]) / sum(abs2, direction)
    kz2 * cos2 + (1 - kz2) / 2 * (1 - cos2)
end

# The recoil second moment of one decay operator per the `recoil_moment`
# setting: exact from the radiation pattern of its emitted component, the
# isotropic 1/3, or a constant.
function decay_recoil_moment(setting::Symbol, label::DecayLabel, mode::MotionalMode)
    if setting === :exact
        recoil_moment(multipole_rank(label.lower, label.upper), label.q, mode.direction)
    elseif setting === :isotropic
        recoil_moment(0, 0, mode.direction)
    else
        throw(
            ArgumentError(
                "recoil_moment must be :exact, :isotropic or a number, got $setting",
            ),
        )
    end
end
decay_recoil_moment(setting::Real, ::DecayLabel, ::MotionalMode) = Float64(setting)
emission_kind(setting::Symbol) = setting
emission_kind(::Real) = :constant

"""
    motional_coupling(species, scheme::LaserScheme, model::LindbladModel, mode; recoil_moment = :exact)
    motional_coupling(species, scheme, model, modes::AbstractVector) -> Vector{MotionalCoupling}

Derives the [`MotionalCoupling`](@ref) of one [`MotionalMode`](@ref) (or of each
of several) to the internal dynamics of `model`, built from `scheme`: the
projected Lamb–Dicke factor of each beam, the sideband Hamiltonian
``\\sum_l η_l \\, i (C_l - C_l^†)/2`` (its beat-note harmonics kept alongside),
and the recoil operators ``\\sqrt{α} η_0 L`` of the spontaneous-emission jump
operators. `recoil_moment = :exact` evaluates ``α`` per emitted component from
the radiation patterns ([`recoil_moment`](@ref)), `:isotropic` takes uniform
emission (``α = 1/3``), and a number uses that constant for every decay (e.g.
the often-quoted ``2/5``) — in which case only this first-order model is
available, cf. [`motional_model`](@ref).
"""
function motional_coupling(
    species,
    scheme::LaserScheme,
    model::LindbladModel,
    mode::MotionalMode;
    recoil_moment=:exact,
)
    basis = scheme.basis
    num_states = length(basis)
    length(model.couplings) == length(scheme.beams) || throw(
        ArgumentError(
            "The model carries $(length(model.couplings)) beam couplings but the " *
            "scheme $(length(scheme.beams)) beams; derive the model from this scheme",
        ),
    )

    # Sideband Hamiltonian: coefficient of (a + a†) from e^{ik·r} ≈ 1 + iη(a + a†),
    # i.e. iη C/2 + h.c. per beam, split into static and beat-note parts like the
    # carrier couplings of the model.
    η = [lamb_dicke(species, mode, beam) for beam in scheme.beams]
    H_sb = zeros(eltype(model.hamiltonian), num_states, num_states)
    harmonics = Tuple{typeof(1.0 * ANGULAR_UNIT),typeof(H_sb)}[]
    for (b, C) in enumerate(model.couplings)
        static, parts = split_beam_coupling(model.frame, b, C)
        H_sb .+= (η[b] * im) .* (static .- static') ./ 2
        for (w, part) in parts
            X = (η[b] * im) .* part ./ 2
            w > zero(w) ? push!(harmonics, (w, X)) : push!(harmonics, (-w, Matrix(X')))
        end
    end
    harmonics = merge_harmonics(harmonics)

    # Recoil: the unprojected Lamb–Dicke parameter of each decay transition,
    # weighted by the emission-pattern second moment along the mode.
    recoil_indices = findall(l -> l isa DecayLabel, model.jump_labels)
    η0 = Float64[]
    α = Float64[]
    ops = similar(model.jump_operators, 0)
    for k in recoil_indices
        label = model.jump_labels[k]
        k_photon = transition_frequency(species, label.lower, label.upper) / u"c"
        η0_k = lamb_dicke(species, mode, k_photon .* mode.direction; projected=false)
        α_k = decay_recoil_moment(recoil_moment, label, mode)
        push!(η0, η0_k)
        push!(α, α_k)
        push!(ops, (sqrt(α_k) * η0_k) .* model.jump_operators[k])
    end

    MotionalCoupling(
        mode,
        η,
        η0,
        α,
        emission_kind(recoil_moment),
        H_sb,
        harmonics,
        ops,
        recoil_indices,
    )
end

motional_coupling(
    species,
    scheme::LaserScheme,
    model::LindbladModel,
    modes::AbstractVector;
    recoil_moment=:exact,
) = [motional_coupling(species, scheme, model, mode; recoil_moment) for mode in modes]


# --- Beyond first order in the Lamb–Dicke parameters ------------------------------

"""
Validates a Lamb–Dicke expansion order: a positive integer, or `Inf` for the
exact operators.
"""
function validate_lamb_dicke_order(order)
    (order isa Real && order >= 1 && (isinf(order) || isinteger(order))) || throw(
        ArgumentError(
            "The Lamb–Dicke order must be a positive integer or Inf, got $order",
        ),
    )
    nothing
end

"""
    displacement_elements(η, num_fock; order = Inf) -> Matrix{ComplexF64}

Returns the matrix elements ``⟨n| D(η) |m⟩``, ``n, m = 0 … \\mathtt{num\\_fock} - 1``,
of the displacement operator ``D(η) = e^{iη(a + a^†)}`` — the motional factor of
the plane wave ``e^{ikx}`` for a mode with Lamb–Dicke factor ``η``, cf.
[`lamb_dicke`](@ref) — or, for a finite `order` ``k``, of its Taylor
truncation ``\\sum_{j ≤ k} (iη(a + a^†))^j / j!`` (`order = 1`: the usual
Lamb–Dicke-regime ``1 + iη(a + a^†)``). The elements are those of the
untruncated operators, evaluated exactly rather than in the truncated Fock
space: in closed form through the generalised Laguerre polynomials,

```math
⟨n| D(η) |m⟩ = e^{-η^2/2} (iη)^{|n - m|} \\sqrt{n_<! / n_>!} \\,
L_{n_<}^{(|n - m|)}(η^2)
```

(Wineland et al., J. Res. Natl. Inst. Stand. Technol. **103**, 259 (1998);
Kulosa et al., New J. Phys. **25**, 053008 (2023), Eq. (7)), evaluated through
the three-term recurrence of the Laguerre polynomials; a finite-order series
is formed in an enlarged Fock space. The red-sideband elements ``⟨n-1|D|n⟩``
saturate and change sign at ``n ≈ 3.67/η^2``, which no finite order captures.
"""
function displacement_elements(η::Real, num_fock::Integer; order::Real=Inf)
    validate_lamb_dicke_order(order)
    num_fock >= 1 || throw(ArgumentError("num_fock must be positive, got $num_fock"))
    if isinf(order)
        C = zeros(ComplexF64, num_fock, num_fock)
        x = η^2
        for d in 0:(num_fock-1)
            # L_k^{(d)}(x) for k = 0 … num_fock − 1 − d (k = n<, n> = k + d)
            L_prev, L = 0.0, 1.0
            for k in 0:(num_fock-1-d)
                if k == 1
                    L_prev, L = L, 1 + d - x
                elseif k > 1
                    L_prev, L = L, ((2k - 1 + d - x) * L - (k - 1 + d) * L_prev) / k
                end
                # √(k!/(k+d)!) η^d, in logs to stay finite for large d
                mag =
                    d == 0 ? 1.0 :
                    exp(d * log(abs(η)) - 0.5 * sum(log, (k+1):(k+d); init=0.0))
                value = exp(-x / 2) * mag * L * (im * sign(η))^d
                C[k+1+d, k+1] = value
                C[k+1, k+d+1] = value
            end
        end
        return C
    end
    Matrix(joint_displacement([η], [num_fock]; order))
end

"""
Returns the (sparse) matrix of the joint displacement operator
``e^{i \\sum_m η_m (a_m + a_m^†)}`` of several modes over their Fock spaces
(last mode fastest), or of its Taylor truncation to the given total `order`
(exact elements, evaluated in an enlarged space). Elements below `cutoff`
times the largest are dropped.
"""
function joint_displacement(ηs, num_fock; order::Real=Inf, cutoff::Real=1e-14)
    validate_lamb_dicke_order(order)
    length(ηs) == length(num_fock) ||
        throw(ArgumentError("One Lamb–Dicke factor per mode is required"))
    D = if isinf(order)
        factors =
            [sparse(displacement_elements(η, nf)) for (η, nf) in zip(ηs, num_fock)]
        reduce(kron, factors; init=sparse(ComplexF64(1) * I, 1, 1))
    else
        k = Int(order)
        enlarged = [nf + k for nf in num_fock]
        # X = Σ_m η_m (a_m + a_m†) on the enlarged joint space
        X = spzeros(ComplexF64, prod(enlarged), prod(enlarged))
        for (m, (η, ne)) in enumerate(zip(ηs, enlarged))
            x = spdiagm(1 => sqrt.(1.0:(ne-1)), -1 => sqrt.(1.0:(ne-1)))
            factors = Any[sparse(1.0I, n, n) for n in enlarged]
            factors[m] = x
            X += η .* reduce(kron, factors)
        end
        term = sparse(ComplexF64(1) * I, size(X)...)
        total = copy(term)
        for j in 1:k
            term = (term * X) .* (im / j)
            total += term
        end
        # restrict to the original Fock spaces
        keep = [
            r for r in 1:prod(enlarged) if all(fock_numbers(r, enlarged) .< num_fock)
        ]
        total[keep, keep]
    end
    droptol!(D, cutoff * maximum(abs, D; init=0.0))
    D
end

# Fock numbers (0-based) of the 1-based full index `r` of a joint space with
# `num_fock[m]` states per mode, last mode fastest (the `tensor` convention of
# the solvers; a preceding internal index, if any, is the slowest).
function fock_numbers(r::Integer, num_fock)
    ns = zeros(Int, length(num_fock))
    x = r - 1
    for m in length(num_fock):-1:1
        ns[m] = x % num_fock[m]
        x ÷= num_fock[m]
    end
    ns
end

# Emitted intensity of the component (rank, q) into the direction with polar
# cosine c (quantisation axis ẑ), up to a constant.
function pattern_value(rank, q, c)
    coefficients = radiation_pattern(rank, q)
    sum(a * c^(i - 1) for (i, a) in enumerate(coefficients))
end

"""
    emission_rule(rank, q, direction; nodes = 12) -> (h, w)
    emission_rule(rank, q, directions::AbstractVector{<:AbstractVector}; nodes = 6) -> (h, w)

Returns a quadrature rule — nodes `h` and weights `w` with ``\\sum_j w_j = 1`` —
for the distribution of the projection ``h = \\hat{k} ⋅ \\vec{e}_m`` of the
direction of a photon emitted in the spherical component `q` of an electric
multipole transition of the given `rank` (`rank = 0`: isotropic emission) onto
the unit mode `direction`, quantisation axis along ẑ:
``∫ K(h) f(h) \\, dh ≈ \\sum_j w_j f(h_j)``. Its second moment is
[`recoil_moment`](@ref); the rule is what the recoil kicks ``e^{-i η_0 h (a +
a^†)}`` of [`motional_model`](@ref) are averaged with beyond first order in
the Lamb–Dicke parameter (Kulosa et al., New J. Phys. **25**, 053008 (2023),
Eq. (A.3)).

For one direction the result is the Gauss rule of the measure ``K(h) \\, dh``
with `nodes` points, exact for polynomials of degree ``2 \\, \\mathtt{nodes} -
1``, obtained by the Stieltjes procedure (Lanczos with full
reorthogonalisation) and the Golub–Welsch eigenvalue method on a fine product
quadrature of the sphere (Gautschi, *Orthogonal Polynomials: Computation and
Approximation*, Oxford (2004), Secs. 2.2 and 3.1). For several `directions`
the nodes are the joint projections of a product rule on the sphere — row
``j`` of the matrix `h` holds ``\\hat{k}_j ⋅ \\vec{e}_m`` for every direction
``m`` — with `nodes` Gauss–Legendre points in ``\\cos θ`` and twice as many
equidistant azimuths, exact for integrands polynomial of degree up to
``2 \\, \\mathtt{nodes} - 5`` in the direction cosines (the pattern itself
being of degree four).
"""
function emission_rule(
    rank::Integer,
    q::Integer,
    direction::AbstractVector{<:Real};
    nodes::Integer=12,
)
    abs(q) <= rank || throw(ArgumentError("|q| must not exceed the rank"))
    nodes >= 1 || throw(ArgumentError("At least one node is required"))
    if rank == 0
        # K(h) = 1/2 on [−1, 1]: the Gauss–Legendre rule itself.
        z, w = gausslegendre(nodes)
        return z, w ./ 2
    end
    e = direction ./ sqrt(sum(abs2, direction))
    hs, ws = sphere_quadrature(rank, q, [e], 48)
    discrete_gauss_rule(vec(hs), ws, nodes)
end

function emission_rule(
    rank::Integer,
    q::Integer,
    directions::AbstractVector{<:AbstractVector};
    nodes::Integer=6,
)
    abs(q) <= rank || throw(ArgumentError("|q| must not exceed the rank"))
    nodes >= 1 || throw(ArgumentError("At least one node is required"))
    es = [d ./ sqrt(sum(abs2, d)) for d in directions]
    sphere_quadrature(rank, q, es, nodes)
end

# Product quadrature of the sphere weighted by the pattern (rank, q):
# Gauss–Legendre in cos θ (n_z points) × 2 n_z equidistant azimuths. Returns the
# projections (node × direction) and the normalised weights.
function sphere_quadrature(rank, q, es, n_z)
    zs, wz = gausslegendre(n_z)
    n_φ = 2n_z
    H = zeros(n_z * n_φ, length(es))
    w = zeros(n_z * n_φ)
    idx = 0
    for (z, wzk) in zip(zs, wz), j in 1:n_φ
        φ = 2π * (j - 1) / n_φ
        s = sqrt(1 - z^2)
        k = (s * cos(φ), s * sin(φ), z)
        idx += 1
        for (m, e) in enumerate(es)
            H[idx, m] = k[1] * e[1] + k[2] * e[2] + k[3] * e[3]
        end
        w[idx] = wzk * pattern_value(rank, q, z)
    end
    w ./= sum(w)
    H, w
end

# Gauss rule with `num_nodes` points for a discrete measure (points x, weights
# w ≥ 0): Stieltjes procedure (Lanczos with full reorthogonalisation) for the
# recurrence coefficients, Golub–Welsch for nodes and weights.
function discrete_gauss_rule(x, w, num_nodes)
    num_points = length(x)
    num_nodes <= num_points ||
        throw(ArgumentError("The rule cannot have more nodes than support points"))
    Q = zeros(num_points, num_nodes)
    α = zeros(num_nodes)
    β = zeros(num_nodes)
    q = sqrt.(w)
    for k in 1:num_nodes
        Q[:, k] = q
        v = x .* q
        α[k] = dot(q, v)
        for j in 1:k   # full reorthogonalisation
            v .-= dot(Q[:, j], v) .* Q[:, j]
        end
        k == num_nodes && break
        β[k] = norm(v)
        q = v ./ β[k]
    end
    J = SymTridiagonal(α, β[1:(num_nodes-1)])
    E = eigen(J)
    E.values, abs2.(E.vectors[1, :])
end

"""
    fock_truncation(nbar; tail = 1e-4, minimum = 2)

Returns the number of Fock states needed to hold a thermal state of mean
occupation `nbar` up to a population `tail` beyond the truncation: the smallest
``N`` with ``\\sum_{n ≥ N} p_n = (\\bar n / (\\bar n + 1))^N ≤ \\mathtt{tail}``,
and at least `minimum`.
"""
function fock_truncation(nbar::Real; tail::Real=1e-4, minimum::Integer=2)
    nbar >= 0 || throw(ArgumentError("nbar must be non-negative"))
    0 < tail < 1 || throw(ArgumentError("tail must lie in (0, 1)"))
    nbar == 0 && return Int(minimum)
    max(Int(minimum), ceil(Int, log(tail) / log(nbar / (nbar + 1))))
end

export MotionalCoupling, lamb_dicke, recoil_moment, motional_coupling
export displacement_elements, emission_rule, fock_truncation
