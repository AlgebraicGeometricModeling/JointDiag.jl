using Test
using JointDiag, LinearAlgebra

N = 4
E = randn(N,N)

eps = 1.e-4

M1 = E*diagm([1,1,2,2])*inv(E) + eps*randn(N,N)
M2 = E*diagm([-1,1,0,3])*inv(E) + eps*randn(N,N)
M3 = E*diagm([1,-1,-2,3])*inv(E) + eps*randn(N,N)

M = [M1,M2,M3]

Xi, E, Slv = joint_diag(M, NewtonJointDiag())

R0 = [E*diagm(Xi[i,:])*inv(E) for i in 1:length(M)]
@test all(norm.(M-R0) .< 0.1)
