using Test
using JointDiag, LinearAlgebra

N = 4
E0 = randn(N+1,N)
F0 = randn(N,N+2)

M00 = E0*F0
M01 = E0*diagm([1,1,2,2])*F0
M02 = E0*diagm([-1,-1,0,3])*F0

M0 = [M00, M01, M02]

X0, U0, V0, Slv = joint_reduce(M0, JRS(1.e-6, NewtonJointDiag()))

R0 = [U0*diagm(X0[i,:])*V0 for i in 1:size(X0,1)]
@test all(norm.(R0-M0) .< 1e-10)
