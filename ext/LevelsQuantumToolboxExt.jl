# Conversions from Levels.jl objects to QuantumToolbox.jl `QuantumObject`s, so
# that the QuantumToolbox solvers (`mesolve`, `steadystate_fourier`, …) can be
# applied to models built from Levels data.
#
# QuantumToolbox works on dimensionless matrices, so the conversions fix a
# reference time unit (the `time_unit` keyword; µs by default, matching the
# µs⁻¹ normalisation of the dynamics layer): a unitful quantity of pure time
# dimension 𝐓ᵖ is stripped in units of `time_unit`ᵖ (µs⁻¹ for Hamiltonians in
# angular units, µs^(-1/2) for Lindblad jump operators, …), and times passed
# to the solvers are then in `time_unit`.

module LevelsQuantumToolboxExt

using LinearAlgebra: Diagonal
using Unitful
using Levels: StateBasis, stateindex
using Levels.PeriodicDynamics: PeriodicDynamics, DrivenTransition
using QuantumToolbox: QuantumToolbox, Ket, Operator, QuantumObject

# Returns the time power p of a quantity of pure time dimension 𝐓ᵖ (0 for a
# dimensionless one), used to fix the stripping unit `time_unit`ᵖ.
function time_power(x)
    dims = typeof(dimension(x)).parameters[1]
    if isempty(dims)
        return 0
    end
    if !(length(dims) == 1 && dims[1] isa Unitful.Dimension{:Time})
        throw(
            ArgumentError(
                "Only quantities of pure time dimension (rates, angular " *
                "frequencies, √rate jump operators, …) map onto the " *
                "dimensionless time-base convention; got $(dimension(x))",
            ),
        )
    end
    dims[1].power
end

# Validates that the reference time unit is indeed one (e.g. `u"µs"`).
function check_time_unit(time_unit)
    if dimension(time_unit) != Unitful.𝐓
        throw(
            ArgumentError(
                "The reference time unit must be a plain unit " *
                "of time; got $time_unit",
            ),
        )
    end
end

strip_time_units(A::AbstractArray{<:Number}, _) = A
strip_time_units(A::AbstractArray{<:Quantity}, time_unit) =
    ustrip.(time_unit^time_power(first(A)), A)

"""
    QuantumObject(A::AbstractMatrix, basis::StateBasis; time_unit = u"µs")
    QuantumObject(v::AbstractVector, basis::StateBasis; time_unit = u"µs")

Builds the `QuantumToolbox.QuantumObject` operator (or ket) representing the
given matrix (or coefficient vector) over the states of `basis`, in
[`Levels.stateindex`](@ref) order.

Unitful elements of pure time dimension ``𝐓^p`` are stripped in units of
``\\mathtt{time\\_unit}^p`` — with the µs default, µs⁻¹ for Hamiltonians in
the angular convention and µs^(-1/2) for Lindblad jump operators — so times
passed to the QuantumToolbox solvers are in `time_unit`.
"""
function QuantumToolbox.QuantumObject(
    A::AbstractMatrix,
    basis::StateBasis;
    time_unit::Unitful.Units=Unitful.µs,
)
    check_time_unit(time_unit)
    if size(A) != (length(basis), length(basis))
        throw(
            ArgumentError(
                "Matrix size $(size(A)) does not match the basis " *
                "($(length(basis)) states)",
            ),
        )
    end
    QuantumObject(strip_time_units(A, time_unit), Operator(), length(basis))
end

function QuantumToolbox.QuantumObject(
    v::AbstractVector,
    basis::StateBasis;
    time_unit::Unitful.Units=Unitful.µs,
)
    check_time_unit(time_unit)
    if length(v) != length(basis)
        throw(
            ArgumentError(
                "Vector length $(length(v)) does not match the basis " *
                "($(length(basis)) states)",
            ),
        )
    end
    QuantumObject(strip_time_units(v, time_unit), Ket(), length(basis))
end

"""
    basis(b::StateBasis, state)
    basis(b::StateBasis, level, m)

Returns the ket `QuantumObject` for the given state of the [`StateBasis`](@ref)
(cf. [`Levels.stateindex`](@ref) for the accepted state forms).
"""
QuantumToolbox.fock(b::StateBasis, state) =
    QuantumToolbox.fock(length(b), stateindex(b, state) - 1)
QuantumToolbox.fock(b::StateBasis, level, m) =
    QuantumToolbox.fock(length(b), stateindex(b, level, m) - 1)

"""
    projection(b::StateBasis, ket[, bra])

Returns the projection operator ``|ket⟩⟨bra|`` between states of the
[`StateBasis`](@ref) (cf. [`Levels.stateindex`](@ref) for the accepted state
forms); `bra` defaults to `ket`, giving the population projector.
"""
QuantumToolbox.projection(b::StateBasis, ket, bra) =
    QuantumToolbox.projection(length(b), stateindex(b, ket) - 1, stateindex(b, bra) - 1)
QuantumToolbox.projection(b::StateBasis, ket) = QuantumToolbox.projection(b, ket, ket)

function PeriodicDynamics.fourier_hamiltonians(
    dt::DrivenTransition;
    δ=zero(dt.drive_frequency),
    time_unit::Unitful.Units=Unitful.µs,
)
    H_0 = (dt.coupling .+ dt.coupling') ./ 2 .+ Diagonal(dt.frame)
    for i in dt.upper_range
        H_0[i, i] -= δ
    end
    H_p = zero(dt.coupling)
    for drive in dt.drives
        H_p .+= (0.5 * cis(drive.phase)) .* drive.amplitude
    end
    (
        H_0=QuantumObject(H_0, dt.basis; time_unit),
        H_p=QuantumObject(H_p, dt.basis; time_unit),
        H_m=QuantumObject(Matrix(H_p'), dt.basis; time_unit),
        ωd=ustrip(time_unit^-1, dt.drive_frequency),
    )
end

end # module
