export off_matrix, JacobiJointDiag, joint_diag

using LinearAlgebra

mutable struct JacobiJointDiag{T}
    max_iter::Int
    epsilon::T
    info::Dict{Symbol,Any}
end

export getindex, setindex!
function Base.getindex(slv::JacobiJointDiag{T}, s::Symbol) where T
    Base.get(slv.info, s, 0)
end
function Base.setindex!(slv::JacobiJointDiag{T}, v, s::Symbol) where T
    slv.info[s] = v
end

JacobiJointDiag() = JacobiJointDiag(1000, 1.e-10, Dict{Symbol,Any}()) 

# norm of off diagonal terms of a square matrix
function norm_off(M)
    if size(M,1)>1
        return sqrt(sum(abs2(M[i,j]) + abs2(M[j,i]) for i in 1:size(M,1) for j in i+1:size(M,1)))
    else
        return 0.0
    end
end

function off_matrix(M)
    C = copy(M)
    for i in 1:size(M,1) C[i,i]=0 end
    return C
end

function jacobi_transform(D, E, F)
     N = [F*D[i]*E for i in 1:length(D)]
end

function jacobi_gradient(D)

    s = length(D)
    O = off_matrix.(D)

    G = sum(D[i]'*O[i]-O[i]*D[i]' for i in 1:length(D))
end

function off_norm2(D)
    sum(norm_off.(D))
end


function jacobi_descent(D)

    s = length(D)
    O = off_matrix.(D)
    S = sum(O[i]*D[i]'-D[i]'*O[i] for i in 1:s)

    H = copy(S)
    for a in 1:size(S,1)
        for b in 1:size(S,2)
            if a != b
                gamma = sum((D[i][a,a]-D[i][b,b])^2 for i in 1:s)
                #println(a, " ", b," ", gamma,)
                H[a,b] /= gamma
            else
                H[a,b] = 0.0
            end
        end
    end

    return S, H
end

function pencil_line_search!(D, S, H)

    #=
    A = [0.,alpha]
    E = I+alpha*S
    F = inv(E)
    N = [F*D[i]*E for i in 1:length(D)]
    GN = -jacobi_gradient(N)
    V = [dot(S,S), dot(S,GN)]
    =#
    
    #=alpha = dot(S,S)/dot(H,S)
    E = I-alpha*S
    F = inv(E)
    global N = [F*D[i]*E for i in 1:length(D)]
    err = sum(norm_off.(N))
    =#
    
    tau = 0.5
    c = 1.e-4
    
    f0 = off_norm2(D)
    df0 = dot(S,jacobi_gradient(D))

    @assert df0<0
    
    alpha = 1.0
    E = I+alpha*S
    N = jacobi_transform(D,E,inv(E))
    falpha = off_norm2(N)

    t = c*df0

    println(">> f0: ", f0,  "  df0: ", df0, "   t: ", t)

    while falpha > f0 + alpha*t && alpha >1.e-8

        alpha *= tau

        E = I + alpha*S
        N = jacobi_transform(D, E, inv(E))
        
        falpha = off_norm2(N)
        
        println("alpha:  ", alpha,  " f:  ", falpha)

    end
#    print("   alpha: ", alpha, "   err: ", falpha)

    # T *= E
    # println(" V: ",V)
    
    N , falpha
end

#----------------------------------------------------------------------
"""
     joint_diag(M::Vector{Matrix{C}}, Solver::JacobiJointDiag)

Compute the joint diagonalization of an array `M` of square matrices `M[1], ..., M[n]` by applying a Jacobi-type iteration to minimize the off-diagonal norm of the matrices `F*M[i]*E` with the constraint `F*E=I`. It outputs

  - `X` the vectors of eigenvalues, which are the columns of X.
  - `E` the common eigenvectors such that `M[i]*E=E*diagm(X[i,:])`

"""
function joint_diag(M::AbstractVector{<:AbstractMatrix{C}},
                              Slv::JacobiJointDiag) where C
    n  = length(M)
    r  = size(M[1],1)

    N   = Slv.max_iter
    eps = Slv.epsilon

    D  = copy(M)

    err = off_norm2(D)

    print("Off0: ", err)

    Slv[:error] = err
    nit = 0
    S0 = zero(M[1])+I
    
    # err=0.0
    if err > eps
        delta = err
        while nit < N && delta > 1.e-4
            err0 = err
            S1, H = jacobi_descent(D)

            #beta = dot(S1, S1-S0)/dot(S0,S0)  #Polak-Ribiere 

            #beta = 0.0
            #S = S1 + beta*S0

            #println("   S ", norm(S1))
            
           
            D, err = pencil_line_search!(D, S1, H) #E, alpha)

            S0 = S1
            nit+=1

            delta = err0-err
            println("\n", nit," Dist: ", err, "   delta: ", delta)
        end
        Slv[:error] = err
    end

    Slv[:nb_iter] = nit
    
    Xi = fill(zero(S0[1,1]),n,r)
    for i in 1:r
    	for j in 1:n
	    Xi[j,i] = D[j][i,i]
	end
    end
    
    return Xi, S0, Slv
end

