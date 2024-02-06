export joint_diag, RandJointDiag

mutable struct RandJointDiag
end


"""
     joint_diag(M::Vector{Matrix{C}}, Solver::RandJointDiag)

Compute the joint diagonalization of an array `M` of square matrices `M[1],...,M[n]` using a random combination of the matrices and its Schur factorization to get the common eigenvectors.
It outputs

  - `X` the vectors of eigenvalues, which are the columns of `X`.
  - `E` the common eigenvectors such that `M[i]*E=E*diagm(X[i,:])`

"""
function joint_diag(M::Vector{Matrix{C}},
                    Solver::RandJointDiag) where C

    M0 = sum(M[i]*randn() for i in 1:length(M))

    E  = schur(M0).vectors

    X = fill(zero(E[1,1]),length(M),size(M0,1))
    for j in 1:length(M)
        Yj = E\(M[j]*E)
        for i in 1:size(M0,1)
            X[j,i]= Yj[i,i] #(Y[:,i]\Yj[:,i])[1] #D[i,i]
        end
    end

    X, E
end
