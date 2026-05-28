using LinearAlgebra

# Manocha, D. & Demmel, J. Algorithms for intersecting parametric and algebraic curves II: multiple intersections
# Graphical Models and Image Processing, Elsevier, 1995, 57, 81-100
function _schur(M::AbstractMatrix{BigFloat})
    sf, Z = _schur(Float64.(M))
    return sf, convert(Matrix{BigFloat}, Z)
end
function _schur(M::AbstractMatrix{<:Real})
    # `M = Z * T * Z'` and `values` gives the eigenvalues
    sf = LinearAlgebra.schur(M)
    return sf, sf.Z
end

"""
    mutable struct SchurJointDiag <: AbstractSolver
        factorization::Union{Nothing,LinearAlgebra.Schur}
    end


Simultaneous diagonalization of commuting matrices using the method of [CGT97].

[CGT97] Corless, R. M.; Gianni, P. M. & Trager, B. M.
*A reordered Schur factorization method for zero-dimensional polynomial systems with multiple roots*
Proceedings of the 1997 international symposium on Symbolic and algebraic computation, 1997, 133-140
"""
mutable struct SchurJointDiag <: AbstractSolver
    factorization::Union{Nothing,LinearAlgebra.Schur}

    function SchurJointDiag()
        return new(nothing)
    end
end

function joint_diag(
    Ms::AbstractVector{<:AbstractMatrix{T}},
    M0::AbstractMatrix,
    solver::SchurJointDiag,
) where {T<:Real}
    sf, Z = _schur(M0)
    solver.factorization = sf
    return solution_matrix(Ms, Z', Z), Z, solver
end

# If i, i+1 are conjugate pair, then they need to be either both in I or both not in I.
# If one of them is in I and the other is not then LAPACK will consider that both of them are in I.
function condition_number(sf::Schur, I)
    n = length(sf.values)
    select = zeros(LinearAlgebra.BlasInt, n)
    for i in I
        select[i] = 1
    end
    try
        return LinearAlgebra.LAPACK.trsen!('E', 'N', select, copy(sf.T), copy(sf.Z))[4]
    catch err
        if err isa LinearAlgebra.LAPACKException
            ε = eps(real(eltype(sf.values)))
            @warn(
                "`LAPACK.trsen!` throwed an exception for `$(I)` so using default tolerance `$ε`"
            )
            # Not sure why this happens, see `test/schur` for an example
            return ε
        else
            rethrow(err)
        end
    end
end

function cluster_arguments(_, s::SchurJointDiag, ɛ)
    # documentation says that the error on the eigenvalues is `ɛ * norm(T) / condition_number`
    nT = norm(s.factorization.T)
    return I -> ɛ * nT / condition_number(s.factorization, I), s.factorization.values
end
