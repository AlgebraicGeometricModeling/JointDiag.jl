export joint_diag, RandJointDiag


mutable struct RandJointDiag
end

function joint_diag(M::Vector{Matrix{C}},
                    Solver::RandJointDiag) where C
    M0 = sum(M[i]*randn() for i in 1:length(M))
    #t0=time()
    I0 = inv(M0)
    #println("... inv   ", time()-t0, "(s)"); t0=time()
    Mg = I0*M[1]

    E  = eigvecs(Mg)
    #println("... eig   ", time()-t0, "(s)"); t0=time()
    Z  = E\I0

    #t0 = time()
    #F = schurfact(Mg)
    #println("... schur ", time()-t0, "(s)"); t0=time()
    # E = F[:vectors]
    # Z = E'

    X = fill(Complex{Float64}(0.0),length(M),size(M0,1))
    for j in 1:length(M)
        Yj = Z*M[j]*E
        # D = Y\Yj
        for i in 1:size(M0,1)
            X[j,i]= Yj[i,i] #(Y[:,i]\Yj[:,i])[1] #D[i,i]
        end
    end
    X, E
end
