using JointDiag, LinearAlgebra

N = 4
E = randn(N,N)

eps = 1.e-6
M1 = E*diagm([1,1,2,2])*inv(E) + eps*randn(N,N)
M2 = E*diagm([-1,1,0,3])*inv(E)

M = [M1,M2]

Xi, E, Slv = joint_diag(M, NewtonJointDiag())

R0 = [E*diagm(Xi[i,:])*inv(E) for i in 1:length(M)]
norm.(M-R0)
