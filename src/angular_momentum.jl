using LinearAlgebra

"""
    jz_matrix(j)

Returns the matrix of the angular momentum operator ``J_z`` for spin `j`, in units
of ``ħ``, in the basis of its eigenstates ordered by increasing projection
``m = -j, …, j``.
"""
jz_matrix(j) = Diagonal([float(m) for m in (-j):j])

"""
    jplus_matrix(j)

Returns the matrix of the raising operator ``J_+`` for spin `j`, in units of ``ħ``,
in the basis of the ``J_z`` eigenstates ordered by increasing ``m``.

The non-zero elements are ``⟨m + 1|J_+|m⟩ = \\sqrt{j (j + 1) - m (m + 1)}``; the
lowering operator ``J_-`` is the adjoint.
"""
function jplus_matrix(j)
    d = Int(2j + 1)
    jp = zeros(d, d)
    for (i, m) in enumerate((-j):j)
        m == j && continue
        jp[i+1, i] = sqrt(float(j * (j + 1) - m * (m + 1)))
    end
    jp
end

"""
    jx_matrix(j)

Returns the matrix of the angular momentum operator ``J_x`` for spin `j`, in units
of ``ħ``, in the basis of the ``J_z`` eigenstates ordered by increasing ``m``.
"""
jx_matrix(j) = (jplus_matrix(j) + jplus_matrix(j)') / 2

"""
    jy_matrix(j)

Returns the matrix of the angular momentum operator ``J_y`` for spin `j`, in units
of ``ħ``, in the basis of the ``J_z`` eigenstates ordered by increasing ``m``.
"""
jy_matrix(j) = (jplus_matrix(j) - jplus_matrix(j)') / (2im)

"""
    spherical_component(v, q)

Returns the spherical component ``q ∈ \\{-1, 0, 1\\}`` of a Cartesian 3-vector:
``v_0 = v_z``, ``v_{±1} = ∓(v_x ± i v_y)/\\sqrt{2}``.
"""
function spherical_component(v, q)
    if q == 0
        complex(float(v[3]))
    elseif q == 1
        -(v[1] + im * v[2]) / sqrt(2)
    elseif q == -1
        (v[1] - im * v[2]) / sqrt(2)
    else
        throw(ArgumentError("Invalid spherical component index q = $q"))
    end
end

public jz_matrix, jplus_matrix, jx_matrix, jy_matrix, spherical_component
