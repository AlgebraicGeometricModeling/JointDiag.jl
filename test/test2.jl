using JointDiag, LinearAlgebra

N = 4
E0 = randn(N,N)

M1 = E0*diagm([1,1,2,2])*inv(E0)
M2 = E0*diagm([-1,-1,0,3])*inv(E0)

M = [M1,M2]

Xi, E = joint_diag(M, RandJointDiag())

