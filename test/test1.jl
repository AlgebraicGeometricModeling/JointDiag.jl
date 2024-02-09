using Test, JointDiag, LinearAlgebra

N = 4
E = randn(N,N)

X = [ 1  1  2 2 ;
     -1  1  0 3 ;
      1 -1 -2 3 ]

M1 = E*diagm(X[1,:])*inv(E) 
M2 = E*diagm(X[2,:])*inv(E) 
M3 = E*diagm(X[3,:])*inv(E) 

M = [M1,M2,M3]

Xi, E0, Slv = joint_diag(M, RandJointDiag())

@test all(norm.(M-[E0*diagm(Xi[i,:])*inv(E0) for i in 1:length(M)]) .< 0.1)

eps = 1.e-5
M1 = M + eps*[randn(N,N) for i in 1:length(M)]

Xi1, E1, Slv1 = joint_diag(M1, NewtonJointDiag())

@test all(norm.(M-[E1*diagm(Xi1[i,:])*inv(E1) for i in 1:length(M)]) .< 0.1)

