using LinearAlgebra

# If i, i+1 are conjugate pair, then they need to be either both in I or both not in I.
# If one of them is in I and the other is not then LAPACK will consider that both of them are in I.
function condition_number(sf::Schur, I)
    n = length(sf.values)
    select = zeros(LinearAlgebra.BlasInt, n)
    for i in I
        select[i] = 1
    end
    try
        return LinearAlgebra.LAPACK.trsen!(
            'E',
            'N',
            select,
            copy(sf.T),
            copy(sf.Z),
        )[4]
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

# Manocha, D. & Demmel, J. Algorithms for intersecting parametric and algebraic curves II: multiple intersections
# Graphical Models and Image Processing, Elsevier, 1995, 57, 81-100
function clusterordschur(M::AbstractMatrix{<:Real}, ɛ)
    if isempty(M)
        # See bug JuliaLang/julia#...
        return Matrix{float(eltype(M))}(undef, 0, 0), Vector{Int}[]
    else
        return _clusterordschur(M, ɛ)
    end
end
function _clusterordschur(M::AbstractMatrix{BigFloat}, ɛ)
    Z, clusters = _clusterordschur(Float64.(M), ɛ)
    return convert(Matrix{BigFloat}, Z), clusters
end
function _clusterordschur(M::AbstractMatrix{<:Real}, ɛ)
    # M = Z * T * Z' and "values" gives the eigenvalues
    sf = LinearAlgebra.schur(M)
    # documentation says that the error on the eigenvalues is ɛ * norm(T) / condition_number
    nT = norm(sf.T)
    clusters = cluster_eigenvalues(sf.values) do I
        return ɛ * nT / condition_number(sf, I)
    end
    return sf.Z, clusters
end

"""
    struct ReorderedSchurSolver{T,RNGT<:Random.AbstractRNG} <: AbstractSolver
        ɛ::T
        rng::RNGT
    end

Simultaneous diagonalization of commuting matrices using the method of [CGT97].

[CGT97] Corless, R. M.; Gianni, P. M. & Trager, B. M.
*A reordered Schur factorization method for zero-dimensional polynomial systems with multiple roots*
Proceedings of the 1997 international symposium on Symbolic and algebraic computation, 1997, 133-140
"""
struct ReorderedSchurSolver{T} <: AbstractSolver
    ɛ::T
end

function ReorderedSchurSolver{T}() where {T}
    return ReorderedSchurSolver(Base.rtoldefault(real(T)))
end

function joint_diag(
    Ms::AbstractVector{<:AbstractMatrix{T}},
    M0::AbstractMatrix,
    solver::ReorderedSchurSolver,
) where {T<:Real}
    n = length(Ms)
    Z, clusters = clusterordschur(M0, solver.ɛ)
    r = length(clusters)
    X = zeros(T, n, r)
    for k in 1:r
        nk = length(clusters[k])
        for j in clusters[k]
            q = Z[:, j]
            for i in 1:n
                X[i, k] += dot(q, Ms[i] * q) / nk
            end
        end
    end
    return Z, X
end
