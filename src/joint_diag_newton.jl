export norm_off, NewtonJointDiag, joint_diag


using LinearAlgebra

mutable struct NewtonJointDiag{T} <: AbstractSolver
    max_iter::Int
    epsilon::T
    info::Dict{Symbol,Any}
end

export getindex, setindex!
function Base.getindex(slv::NewtonJointDiag{T}, s::Symbol) where T
    Base.get(slv.info, s, 0)
end
function Base.setindex!(slv::NewtonJointDiag{T}, v, s::Symbol) where T
    slv.info[s] = v
end

NewtonJointDiag() = NewtonJointDiag(10, 1.e-10, Dict{Symbol,Any}()) 

# norm of off diagonal terms of a square matrix
function norm_off(M)
    if size(M,1)>1
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

#----------------------------------------------------------------------
"""
     joint_diag(M::Vector{Matrix{C}}, Solver::NewtonJointDiag)

Compute the joint diagonalization of an array `M` of square matrices `M[1], ..., M[n]` by applying a Newton-type iteration to minimize the off-diagonal norm of the matrices `F*M[i]*E` with the constraint `F*E=I`. It outputs

  - `X` the vectors of eigenvalues, which are the columns of X.
  - `E` the common eigenvectors such that `M[i]*E=E*diagm(X[i,:])`

It implements the method described in

[KMY22] Khouja, Rima, Mourrain, Bernard, and Yakoubsohn, Jean-Claude.
*Newton-type methods for simultaneous matrix diagonalization.*
Calcolo 59.4 (2022): 38. doi:10.1007/s10092-022-00484-3, https://hal.science/hal-03390265.

"""
function joint_diag(
    M::AbstractVector{<:AbstractMatrix{C}},
    M0::AbstractMatrix,
    Slv::NewtonJointDiag,
) where C
    n  = length(M)
    r  = size(M[1],1)

    N   = Slv.max_iter
    eps = Slv.epsilon

    E  = eigvecs(M0)

    F  = inv(E)
    
    D  = vcat([Matrix{C}(I,r,r)],[F*M[i]*E for i in eachindex(M)])
    err = sum(norm_off.(D))
    delta = sum(norm.(D))

    #println("Off0: ", err, "    delta: ", err/delta)

    Slv[:error] = err
    nit = 0

    if err/delta > eps
        delta = err
        while nit < N && delta > eps
            err0 = err
            X,Y = newton_joint_diag_iter(D)
            D = [Y*D[i]*X for i in eachindex(D)]
            E = E*X
            F = Y*F
            nit+=1
            err = sum(norm_off.(D))
            delta = err0-err
            #println("Off", nit,": ", err, "   delta: ", delta)
        end
        Slv[:error] = err
    end

    Slv[:nb_iter] = nit
    
    Xi = fill(zero(E[1,1]),n,r)
    for i in 1:r
    	for j in 1:n
	    Xi[j,i] = D[j+1][i,i]/D[1][i,i]
	end
    end
    return Xi, E, Slv
end

