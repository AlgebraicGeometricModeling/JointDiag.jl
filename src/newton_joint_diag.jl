export norm_off, NewtonJointDiag, joint_diag

using LinearAlgebra

mutable struct NewtonJointDiag{T}
    max_iter::Int
    epsilon::T
    error::T
    nb_iter::Int
end

NewtonJointDiag() = NewtonJointDiag(10, 1.e-3, 0.0, 0) 

# norm of off diagonal terms of a square matrix
function norm_off(M)
    if size(M[1],1)>1
        return sqrt(sum(abs2(M[i,j]) + abs2(M[j,i]) for i in 1:size(M,1) for j in i+1:size(M,1)))
    else
        return 0.0
    end
end

function newton_joint_diag_iter(D)
    n = size(D[1],1)
    s = length(D)
    
    X = fill(zero(D[1][1,1]),n,n)
    Y = fill(zero(D[1][1,1]),n,n)

    A = fill(zero(D[1][1,1]),s,2)
    b = fill(zero(D[1][1,1]),s)
    for i in 1:n
        for j in 1:n
            if i != j
                for k in 1:s
                    A[k,1] = D[k][i,i]
                    A[k,2] = D[k][j,j]
                    b[k]   = -D[k][i,j]
                end
                v = A\b
                X[i,j] =  v[1]
                Y[i,j] =  v[2]
            end
        end
    end
    for i in 1:n
        X[i,i]=1
        Y[i,i]=1
    end
    return X, Y
end

function joint_diag(M::Vector{Matrix{C}},
                    Solver::NewtonJointDiag) where C
    n  = length(M)
    r  = size(M[1],1)

    N   = Solver.max_iter
    eps = Solver.epsilon

    M1 = sum(M[i]*randn(Float64) for i in 1:n)
    E  = eigvecs(M1)

    F  = inv(E)
    
    D  = vcat([Matrix{C}(I,r,r)],[F*M[i]*E for i in 1:length(M)])
    err = sum(norm_off.(D))
    delta = sum(norm.(D))
    #println("diag off: ", err)

    Solver.error = err
    nit = 0

    if err/delta > 5.e-2
        delta = err
        while nit < N && delta > eps
            err0 = err
            X,Y = newton_joint_diag_iter(D)
            D = [Y*D[i]*X for i in 1:length(D)]
            E = E*X
            F = Y*F
            nit+=1
            err = sum(norm_off.(D))
            delta = err0-err
            #println("Off", nit,": ", err, "   delta: ", delta)
        end
        Solver.error = err
    end

    Solver.nb_iter = nit
    
    Xi = fill(zero(E[1,1]),n,r)
    for i in 1:r
    	for j in 1:n
	    Xi[j,i] = D[j+1][i,i]/D[1][i,i]
	end
    end
    return Xi, E, Solver
end

