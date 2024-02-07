using Test
using JointDiag, LinearAlgebra

N = 4
E0 = randn(N,N)

M01 = E0*diagm([1,1,2,2])*inv(E0)
M02 = E0*diagm([-1,-1,0,3])*inv(E0)

M0 = [M01,M02]

X0, E0 = joint_diag(M0, RandJointDiag())

R0 = [E0*diagm(X0[i,:])*inv(E0) for i in 1:length(M0)]
@test norm.(M0-R0) ≈ [8.74, 31.36] rtol=1e-3

M11 = [0 0 0 0;
       1 0 0 0;
       0 0 0 0;
       0 0 1 0]

M12 = [0 0 2 0;
       0 0 0 2;
       1 0 0 0;
       0 1 0 0]

M1 = [M11, M12]

X1, E1 = joint_diag(M1, RandJointDiag())

R1 = [E1*diagm(X1[i,:])*inv(E1) for i in 1:length(M1)]
@test norm.(M1-R1) ≈ [√2, √2]