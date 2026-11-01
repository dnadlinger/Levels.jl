using OrdinaryDiffEqVerner: Vern8

"A numerical method for the sideband response of a driven transition."
abstract type SidebandMethod end

"""
    DressedStates(; harmonic_cutoff=8, num_time_samples=256)

Dresses each manifold with a Floquet truncation of harmonics `-K:K`, then
computes the laser coupling to first order on `num_time_samples` times per period.
The cutoff and time sampling should be increased to check convergence.
"""
struct DressedStates <: SidebandMethod
    harmonic_cutoff::Int
    num_time_samples::Int
    function DressedStates(; harmonic_cutoff::Integer=8, num_time_samples::Integer=256)
        harmonic_cutoff >= 0 ||
            throw(ArgumentError("harmonic_cutoff must be non-negative"))
        num_time_samples >= 1 ||
            throw(ArgumentError("num_time_samples must be positive"))
        new(harmonic_cutoff, num_time_samples)
    end
end

"""
    Monodromy(; bracket_scale=8.0, estimate=DressedStates(), alg=Vern8(),
               reltol=1e-10, abstol=1e-12)

Finds the sideband resonance through one-period TDSE propagation and an
avoided-crossing search. Includes the laser coupling nonperturbatively within
the model. `estimate` supplies the search centre, `alg` the ODE integrator.
Numerical accuracy remains controlled by the propagation tolerances and search.
"""
struct Monodromy{A} <: SidebandMethod
    bracket_scale::Float64
    estimate::DressedStates
    alg::A
    reltol::Float64
    abstol::Float64
    function Monodromy(;
        bracket_scale::Real=8.0,
        estimate=DressedStates(),
        alg=Vern8(),
        reltol::Real=1e-10,
        abstol::Real=1e-12,
    )
        bracket_scale > 0 || throw(ArgumentError("bracket_scale must be positive"))
        (reltol > 0 && abstol > 0) ||
            throw(ArgumentError("Tolerances must be positive"))
        new{typeof(alg)}(bracket_scale, estimate, alg, reltol, abstol)
    end
end

"A numerical method for the quasienergy shift under harmonic driving."
abstract type QuasienergyShiftMethod end

"""
    PerturbationTheory(; warn_ratio=20)

Evaluates the second-order Floquet perturbation sum. Warns when a denominator
is smaller than `warn_ratio` times its coupling.
"""
Base.@kwdef struct PerturbationTheory <: QuasienergyShiftMethod
    warn_ratio::Float64 = 20.0
end

"""
    Floquet(; harmonic_cutoff=8, min_overlap=0.5)

Evaluates the nonperturbative quasienergy shift from Floquet dressing on
`-harmonic_cutoff:harmonic_cutoff`. Warns when identification of the dressed
state by its bare-state overlap becomes ambiguous.
"""
struct Floquet <: QuasienergyShiftMethod
    harmonic_cutoff::Int
    min_overlap::Float64
    function Floquet(; harmonic_cutoff::Integer=8, min_overlap::Real=0.5)
        harmonic_cutoff >= 0 ||
            throw(ArgumentError("harmonic_cutoff must be non-negative"))
        0 <= min_overlap <= 1 || throw(ArgumentError("min_overlap must lie in [0,1]"))
        new(harmonic_cutoff, min_overlap)
    end
end

export SidebandMethod, DressedStates, Monodromy
export QuasienergyShiftMethod, PerturbationTheory, Floquet
