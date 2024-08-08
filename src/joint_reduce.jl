export joint_reduce, JRS

using LinearAlgebra

rkf_eps = eps::Float64 -> function (S)
  i :: Int = 1;
  while i<= length(S) && S[i]/S[1] > eps
    i+= 1;
  end
  i-1;
end

rkf_cst = r::Int64 -> function (S) return r end

"""
Describe the Joint Reduce Solver

     - `rank` is a rank function used to determine the rank from the singular values
     - `diag_solver`is the solver used for the joint diagonalization

"""
mutable struct JRS
    rank::Function
    diag_solver::Any
end

JRS() = JRS(rkf_eps(1.e-6), RandJointDiag())
JRS(r::Int, slvr = RandJointDiag()) = JRS(rkf_cst(r),slvr)
JRS(epsilon::Float64, slvr = RandJointDiag()) = JRS(rkf_eps(epsilon),slvr)

"""
     joint_reduce(H::Vector{Matrix{C}}, Solver = JRS())

Compute the joint diagonal reduction of an array `H` of square matrices `H[1],...,H[n]` using a random combination of the matrices, SVD and rank reduction to obtain vectors of eigenvalues X,  left factor U, right factor V such that `H=[U*diagm(X[i,:])*V for i in 1:length(H)]`
It outputs

  - `X` the eigenvalue matrix, which columns are the vectors of common eigenvalues
  - `U` the left factor, which columns are  orthogonal
  - `V` the right factor, which rows are orthogonal

"""
function joint_reduce(H::Vector{Matrix{C}},
                      Slv::JRS = JRS()) where C

    H0 = sum(H[i]*randn() for i in 1:length(H))

    n = length(H)

    U, S, V = LinearAlgebra.svd(H0)       # H0= U*diag(S)*V'
    r = Slv.rank(S)

    Sr  = S[1:r]
    Sri = LinearAlgebra.diagm(inv.(Sr))

    M = [ Sri*(U[:,1:r]')*H[i]*(V[:,1:r]) for i in 1:length(H) ]

    if r > 1
        Xi, E = JointDiag.joint_diag(M, Slv.diag_solver)
    else
        Xi = zeros(C, n, r)
        for i in 1:n
            Xi[i,1] = M[i][1,1]
        end
        E  = ones(C, 1, 1)
    end

    Uxi = (U[:,1:r].*Sr')*E
    Vxi = (E\ V[:,1:r]')

    return Xi, Uxi, Vxi, Slv
end
