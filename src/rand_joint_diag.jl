export joint_diag, RandJointDiag


mutable struct RandJointDiag
end

function joint_diag(M::Vector{Matrix{C}},
                    Solver::RandJointDiag) where C

    M0 = sum(M[i]*randn() for i in 1:length(M))

    E  = schur(M0).vectors

    X = fill(zero(E[1,1]),length(M),size(M0,1))
    for j in 1:length(M)
        Yj = E\(M[j]*E)
        for i in 1:size(M0,1)
            X[j,i]= Yj[i,i] #(Y[:,i]\Yj[:,i])[1] #D[i,i]
        end
    end

    X, E
end
