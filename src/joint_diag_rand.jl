export joint_diag, RandJointDiag

mutable struct RandJointDiag <: AbstractSolver end

"""
     joint_diag(M::Vector{Matrix{C}}, Solver::RandJointDiag)

Compute the joint diagonalization of an array `M` of square matrices `M[1],...,M[n]` using a random combination of the matrices and its Schur factorization to get the common eigenvectors.
It outputs

  - `X` the vectors of eigenvalues, which are the columns of `X`.
  - `E` the common eigenvectors such that `M[i]*E=E*diagm(X[i,:])`

"""
function joint_diag(M::Vector{Matrix{C}}, M0::AbstractMatrix, Slv::RandJointDiag) where {C}
    E = eigvecs(M0)
    F = inv(E)

    X = fill(zero(E[1, 1]), length(M), size(M0, 1))
    for j in 1:length(M)
        Yj = F * (M[j] * E)
        for i in axes(M0, 1)
            X[j, i] = Yj[i, i] #(Y[:,i]\Yj[:,i])[1] #D[i,i]
        end
    end

    X, E, Slv
end
