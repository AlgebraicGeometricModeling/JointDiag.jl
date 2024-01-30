using JointDiag, LinearAlgebra

N = 4
E = randn(N,N)

M1 = E*diagm([1,1,2,2])*inv(E)
M2 = E*diagm([-1,1,0,3])*inv(E)


Xi, E = joint_diag([M1,M2], NewtonJointDiag())
