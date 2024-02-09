using Test, JointDiag, LinearAlgebra

N = 4
E0 = randn(N,N)

function mmt(X) 
    [sum(X[1,k] for k in 1:N),
     sum(X[2,k] for k in 1:N),
     sum(X[1,k]^2 for k in 1:N),
     sum(X[1,k]*X[2,k]^2 for k in 1:N),
     sum(X[2,k]^2 for k in 1:N),
     sum(X[1,k]^2*X[2,k] for k in 1:N),
     sum(X[1,k]*X[2,k]^2 for k in 1:N),
     ]
end


Xi = [  1  1 2 2 ;
       -1 -1 0 3 ]

M01 = E0*diagm(Xi[1,:])*inv(E0)
M02 = E0*diagm(Xi[2,:])*inv(E0)

M0 = [M01,M02]

X0, E0 = joint_diag(M0, RandJointDiag())

#R0 = [E0*diagm(X0[i,:])*inv(E0) for i in 1:length(M0)]

@test mmt(X0) ≈ mmt(Xi) rtol=1e-3

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

#R1 = [E1*diagm(X1[i,:])*inv(E1) for i in 1:length(M1)]

Xi1 = [ 0   0   0   0 ;
        √2 √2 -√2 -√2 ]
@test mmt(X1) ≈ mmt(Xi1) rtol=1e-3 

