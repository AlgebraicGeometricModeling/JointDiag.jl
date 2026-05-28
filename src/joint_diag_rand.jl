export joint_diag, EigenJointDiag

mutable struct EigenJointDiag <: AbstractSolver end

"""
     joint_diag(M::Vector{Matrix{C}}, Solver::EigenJointDiag)

Compute the joint diagonalization of an array `M` of square matrices `M[1],...,M[n]` using a random combination of the matrices and its Schur factorization to get the common eigenvectors.
It outputs

  - `X` the vectors of eigenvalues, which are the columns of `X`.
  - `E` the common eigenvectors such that `M[i]*E=E*diagm(X[i,:])`

"""
function joint_diag(M::Vector{Matrix{C}}, M0::AbstractMatrix, Slv::EigenJointDiag) where {C}
    E = eigvecs(M0)
    F = inv(E)
    return solution_matrix(M, F, E), E, Slv
end
