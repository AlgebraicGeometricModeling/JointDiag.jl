using Test, JointDiag, LinearAlgebra

includet("../src/joint_diag_jacobi.jl")

function mmt(X, N = size(X,2))
    
    [sum(X[1,k] for k in 1:N),
     sum(X[2,k] for k in 1:N),
     sum(X[1,k]^2 for k in 1:N),
     sum(X[1,k]*X[2,k]^2 for k in 1:N),
     sum(X[2,k]^2 for k in 1:N),
     sum(X[1,k]^2*X[2,k] for k in 1:N),
     sum(X[1,k]*X[2,k]^2 for k in 1:N),
     ]
end


N = 4
E = randn(N,N)

X = [ 1  1  2 2 ;
     -1  1  0 3 ;
      1 -1 -2 3 ]
    

eps = 1.e-3

M = [E*diagm(X[i,:])*inv(E) + eps*randn(N,N) for i in 1:size(X,1)]

Xi1, E1, Slv1 = joint_diag(M, JacobiJointDiag())

res1 = norm.(M-[E1*diagm(Xi1[i,:])*inv(E1) for i in 1:length(M)])
println("Jacobi: ", res1,  "  ", norm(mmt(X) - mmt(Xi1)))

@test norm(mmt(X) - mmt(Xi1)) < 100

Xi, E, Slv = joint_diag(M, NewtonJointDiag())

res0 = norm.(M-[E*diagm(Xi[i,:])*inv(E) for i in 1:length(M)])
println("Newton: ", res0, "  ", norm(mmt(X) - mmt(Xi)))      

@test norm(mmt(X) - mmt(Xi)) < 100
