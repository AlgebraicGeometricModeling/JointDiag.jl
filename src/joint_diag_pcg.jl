# 1. Import libraries
using DynamicPolynomials, TensorDec, Random, LinearAlgebra, MultivariateSeries, Plots
###################################


###################################
# 2. Define function
## 2.1 Main function: jeigen_pcg
function jeigen_pcg(As::Vector{Matrix{Float64}}; maxiter=500, tol=1e-6, regul=1e-2, init="eye", verbose=true)
    """
    Objective: the function computes joint diagonalization B of a collection of same size matrices 
    (Note that: the columns of B are right eigenvectors, and the rows of inv(B) are left eigenvectors of Aq)

    Args:
        As: vector of Q square matrices of size (m,m)
        maxiter: maximum iterations 
        tol: tolerance 
        regul: regularizer
        init: initial choice for B
        verbose: controlling parameter

    Outputs:
        B: matrix of size (m,m) such that each inv(B) * A[:,:,q] * B is as diagonal as possible.
        trajectories: dictionary track information
        B_list: list of B at each iteration

    Method:
        Minimize the off(As)
        1. compute descent direction: R = - G(T)
        2. do backtracking linesearch: alpha = argmin(off(I + x * R + 0.5 * x^2 * R * R) * T )
        3. update T = (I + x * R + 0.5 * x^2 * R * R) * T
    """

    # Step 1. Pre-processing
    As = cat(As..., dims=3) # stack As to get an array of size (m,m,Q)
    tstart = time()

    m, jk, Q = size(As)
    @assert m == jk "square matrices, please"

    if Q == 1
        D, B = eigen(As[:,:,1])
        # B = normalize_columns(B)
        trajectories = Dict{Symbol, Any}()  # Initialize trajectories dictionary
        trajectories[:moveon] = zeros(0)
        trajectories[:moveoff] = zeros(0)
        trajectories[:jdcrit] = zeros(0)
        trajectories[:numbt] = zeros(0)
        trajectories[:gradoff] = zeros(0)
        trajectories[:gradon] = zeros(0)
        return B, trajectories
    end

    regul_on = 1e-6
    inner_balancing = false

    if verbose
        println("joint diagonlization of $Q matrices of size $m")
    end

    if init == "eye"
        B = Matrix{Float64}(I, m, m)
    elseif init == "rand"
        B = randn(Complex{Float64}, m, m)
    elseif init == "algebra"
        # Perform eigen decomposition and initialize B accordingly
        D, B = eigen(As[:,:,1]-As[:,:,2]) #????
        # B = normalize_columns(B)
    else
        error("Unknown init option")
    end

    As = updateA(As, B, inv(B))
    o = off(As)
    # Bup = normalize_columns(inv(B))
    B_list = Any[B]

    if verbose
        println("JD crit. after init : $o")
    end

    # L, As, iterbal = jd_balance_QN(As, false)
    L, As, iterbal = jd_balance_QN(As; verbose= false)
    B = diagm(L) * B
    o = off(As)

    if verbose
        println("JD crit. after balancing powers : $o")
    end

    # Step 2. Main loop
    ## Step 2.1 Initialize parameters
    jdcrit = zeros(maxiter)
    moveoff = zeros(maxiter)
    moveon = zeros(maxiter)
    numbt = zeros(maxiter)
    gradoff = zeros(maxiter)
    gradon = zeros(maxiter)

    G = zeros(m)
    H = zeros(m)
    D = zeros(m)
    Hc = zeros(m)
    T = zeros(m)
    S = zeros(m)
    R = zeros(m)

    cg_reset = m * m
    k = 0

    ## Step 2.2 Compute descent gradient
    G, H, Hc = der12(As, isCalculateHessian=true)  # derivatives
    R = -G
    S = applyprecon(R, H, Hc, regul, regul_on)
    D = S
    delta_new = real(dot(R[:],  D[:]))
    g0 = 2 * real(dot(G[:],  D[:]))

    trajectories = Dict{Symbol, Any}()
    ## Step 2.3 Do backtracking linesearch and update
    let count_iter = 0
        for iter in 1:maxiter
            # Line search
            As, onew, T, nbt, report = linesearch(As, o, D, g0)
            if nbt == -1
                if verbose
                    println("Max backtrack reached.")
                end
                count_iter = iter
                break
            end
            
            # Update
            B = T * B
            gain = o - onew
            o = onew
            # Bup = normalize_columns(inv(B))
            push!(B_list, B)
            
            # New derivatives
            G, H, Hc = der12(As, isCalculateHessian=true)
            gradoff[iter] = norm(G - diagm(diag(G)))
            gradon[iter] = norm(diagm(diag(G)))
            
            # Stopping based on max relative move
            maxmove = maximum(abs.(T - diagm(diag(T))))
            if maxmove < tol
                if verbose
                    println("Tolerance reached: no relative move larger than $tol")
                end
                count_iter = iter
                break
            end
            
            # Preconditioned conjugation with Polak-Ribiere
            R = -G
            delta_old = delta_new
            delta_mid = real(dot(R[:], S[:]))
            S = applyprecon(R, H, Hc, regul, regul_on)
            delta_new = real(dot(R[:], S[:]))
            beta = (delta_new - delta_mid) / delta_old
            D = S + beta * D
            
            g0 = 2 * real(dot(D[:], G[:]))
            
            CG_RESET = (g0 > 0) || (beta < 0) || (k == cg_reset)
            if CG_RESET
                D = S
                k = 0
                g0 = 2 * real(dot(D[:],  G[:]))
                println("CG reset")
            else
                k += 1
            end
            
            if inner_balancing
                L, As, iterbal = jd_balance_QN(As)
                B = diagm(L) * B
            end
            
            jdcrit[iter] = onew
            numbt[iter] = nbt
            moveon[iter] = sqrt(sum((abs.(diag(T)) .- 1).^2))
            moveoff[iter] = sqrt(sum(abs.(T - diagm(diag(T))) .^ 2))
            
            if verbose
                println("$iter | jd $o (- $gain) | on/off = $(moveon[iter]) / $(moveoff[iter]) | grad on/off = $(gradon[iter]) / $(gradoff[iter]) | nbt = $nbt $report | CG: G0 = $g0 beta = $beta $CG_RESET")
            end
            count_iter = iter
        end
    
    if count_iter == maxiter
        if verbose
            println("Max number of iterations reached")
        end
    end
    
    if verbose
        dur = time() - tstart
        println("JD criterion = $o. Relative criterion = $(off_relative(As)) [in $dur secs]")
    end
    
    # Final normalization for one global scale and phase
    B = B / tr(B)
    # Binv = copy(B)
    B = inv(B)
    # B = normalize_columns(inv(B))
    
    F = zeros(Q,m)
    for i=1:Q
        #        F[i,:] = vec(diag(B * As[:,:,i] * inv(B)))
        F[i,:] = diag(As[:,:,i])    
    end


    # Trim trajectories
    keep = 1:(count_iter - 1)
    trajectories[:moveon] = moveon[keep]
    trajectories[:moveoff] = moveoff[keep]
    trajectories[:jdcrit] = jdcrit[keep]
    trajectories[:numbt] = numbt[keep]
    trajectories[:gradoff] = gradoff[keep]
    trajectories[:gradon] = gradon[keep]

    return F, B, trajectories, B_list, o
    end
end
###################################


###################################
## 2.2 Local functions
### 2.2.1 off
function off(As)
    """
    Arg:
        As: an array of size (m,m,Q)

    Output:
        o: the criterion, sum of off-diagonal Frobenius norms of A[:,:,q]
    """

    Q = size(As,3)
    o = 0
    for q=1:Q
        M = As[:, :, q]
        M = M - diagm(diag(M))
        o += dot(M[:], M[:])
    end
    return o
end
###################################


### 2.2.2 off_relative 
#### (be not used in the main function)
function off_relative(As)
    """
    Arg:
        As: an array of size (m,m,Q)

    Output:
        o: normalization of off function by its diagonal
    """
    m, jk, Q = size(As)
    
    # Initialize variables
    D = zeros(m)
    M = zeros(m, m)
    o = 0
    d = 0
    
    # Iterate over each matrix in As
    for q = 1:Q
        M = As[:, :, q]
        D = diag(M)
        M = M - diagm(diag(M))
        o += dot(M[:], M[:])
        d += dot(D, D)
    end
    
    # Calculate normalized version of JD criterion
    or = sqrt(o / d / m / Q)
    return or
end
###################################

### 2.2.3 updateA
function updateA(A, T, iT)
    """
    This function is to update the array As.
    In the context: update inv(T) * A * T

    Args:
        As: array of size (m,m,Q)
        T, iT: matrices of size (m,m)

    Output:
        Anew = T * A * iT component-wise
    """
    Anew = copy(A)
    Q = size(A, 3)
    for q=1:Q
        Anew[:, :, q] = T * A[:, :, q] * iT
    end
    return Anew
end
###################################


### 2.2.4 normalize_columns
function normalize_columns(A::AbstractMatrix)
    for j in axes(A, 2)  # Iterate over columns
        A[:, j] .= A[:, j] ./ A[1, j]  # Element-wise division
    end
    return A
end
###################################


### 2.2.5 der12
function der12(AA; isCalculateHessian::Bool=false)
    # 1st and 2nd-order quantities
    """
    Objective: compute Gradient and Hessian of off(AA)

    Args:
        AA: array of size (m,m,Q)
        isCalculateHessian:
            true: compute Gradient, Hessian and conjugate Hessian
            false: compute only Gradient
    
    Output:
        G: Gradient
        H: Hessian
        Hc: conjugate Hessian
    """

    m, jk, Q = size(AA)
    G = zeros(m, m)
    
    for q = 1:Q
        A = AA[:, :, q]
        dA = diag(A)
        
        Ao = A - diagm(diag(A))
        G += Ao * A' - A' * Ao
    end
    
    if isCalculateHessian
        H = zeros(m, m)
        Hc = zeros(m, m)

        for q = 1:Q
            A = AA[:, :, q]
            dA = diag(A)

            D = (A .* conj(A))'
            r = sum(D, dims=1)
            c = sum(D, dims=2)
            H +=  - 2 * (D - diagm(diag(D)) + real(dA * dA')) .+ r .+ c
            Hc += diagm(vec(r)) + diagm(vec(c)) - D - D'
        end

        return G, H, Hc
    else
        return G
    end
end
###################################



### 2.2.6 applyprecon
function applyprecon(R, H, Hc, regul, regul_on)
    """
    The function is to apply a preconditioner to R = -G (descent gradient)

    Args:
        R = -G : descent gradient
        H: hessian
        Hc: conjugate hessian
        regul: regularizer from main function
        regul_on: active regularizer
    
    Output:
        S: adjusted direction 
    """
    if false
        S = R / maximum(abs(H[:]))
        return S
    end

    m = size(R, 1)
    S = R ./ (H .+ regul)

    if true
        e = (Hc + ones(m, m) / m + regul_on * Matrix{Float64}(I, m, m)) \ diag(R)
    else
        e = diag(R) ./ diag(Hc .+ regul_on )
    end

    S += diagm(e / 2 - diag(S))
    return S
end
###################################



### 2.2.7 linesearch
function linesearch(As, o, D, g0)
    """
    Backtracking linesearch using first Wolfe condition:
    f(x + alpha) =< f(x) + c_wolfe * alpha * <p, G>
    where p descent direction and G is gradient
    Here we use c_wolfe = 0.1

    Args: 
        As: array of size (m,m,Q)
        o = off(As): objective function
        D: adjusted descent gradient (after using applyprecon)
        g0 = 2 * <D, G>

    Outputs:
        Anew: updated As by updated T 
        onew: updated objective function 
        T = (I + x * D + 0.5 * x^2 * D * D) * T
        ibt_count: number of excuted iterations
        report: report each iteration

    Method:
        1. check Wolfe condition
        2. update step size
        3. try secant step, update 
    """
    # line search parameters
    backtrackfactor = 2
    maxbacktrack = 20
    c_wolfe = 0.1

    D2 = D * D
    m = size(D, 1)
    Anew = zeros(size(As))

    # 1: backtracking
    alpha = 1
    btsuccess = false
    ibt_count = 0
    for ibt = 1:maxbacktrack
        T = Matrix{Float64}(I, m, m) + alpha * D + 0.5 * alpha^2 * D2
        Anew = updateA(As, T, inv(T))
        onew = off(Anew)
        goodstep = onew < (o + alpha * c_wolfe * g0)
        # ibt_count +=1
        if goodstep
            btsuccess = true
            ibt -= 1
            ibt_count = ibt
            break
        end
        alpha /= backtrackfactor
        # ibt_count = ibt
    end

    if !btsuccess
        ibt_count = -1  # used as a flag
        report = "Backtrack failed"
        Anew = As
        onew = o
        T = Matrix{Float64}(I, m, m)
        return Anew, onew, T, ibt_count, report
    else
        gain_bt = o - onew
    end

    # 2: then try one secant step, starting at backtrack value
    Gnew = der12(Anew)
    g1 = 2 * real(dot(Gnew[:], D[:]))

    wrtbk = g0 / (g0 - g1) # step wrt backtrack step by the secant rule
    alpha_sec = alpha * wrtbk
    Tsec = Matrix{Float64}(I, m, m) + alpha_sec * D + 0.5 * alpha_sec^2 * D2
    Asec = updateA(As, Tsec, inv(Tsec))
    osec = off(Asec)
    gainsec = o - osec

    # check that secant step actually improves ovr backtrack (much needed during first iterations)
    if gainsec > gain_bt # use secant result
        gain = gainsec
        alpha = alpha_sec
        Anew = Asec
        onew = osec
        T = Tsec
        Gnew = der12(Anew)
        secimprove = gainsec / gain_bt
    else # fall back to backtrack values
        gain = gain_bt
        wrtbk = 1
        secimprove = 0
    end

    # relative gradient decrease during line search
    g2 = 2 * real(dot(Gnew[:], D[:]))

    wc2 = g2 / g0

    report = "|LS: nbt=$(ibt_count) g=$(secimprove) rstep=$(wrtbk) wc2=$(wc2)"
    return Anew, onew, T, ibt_count, report
end
###################################




### 2.2.8 jd_balance_QN
function jd_balance_QN(AA; verbose=false)
    """
    This function is to minimize the JD criterion wrt diagonal transforms, that is, to fix scales, which amounts to balancing upper and lower energy
    
    Args:
        AA: array of size (m,m,Q)
        verbose: controlling

    Outputs:
        L: vector containing the new diagonal values after balancing
        AA: array after balancing the diagonals according to L
        count_iter: number of iterations

    Method:
        Updates the diagonal values in vector L and matrix AA based on the Quasi-Newton method 
    """

    m, jk, Q = size(AA)

    N = ones(m, m) / m

    tol = 1e-8
    itermax = 10

    B = zeros(m, m)
    for q = 1:Q
        Aq = AA[:, :, q]
        B += (Aq .* conj(Aq))
    end

    L = ones(m)
    count_iter = 0
    for iter in 1:itermax
        o1 = sum(sum(B - diagm(diag(B))))

        r = sum(B, dims=1)
        c = sum(B, dims=2)
        g = vec(c) - vec(r)
        H = diagm(vec(r)) + diagm(vec(c)) - B - B'
        d = - (H + N) \ g

        d[isnan.(d)] .= 0 # replace NaN by zero

        d = exp.(d / 2)
        L = L .* d
        d = d .* d  # square to act on B
        B = d .* B .* (1.0 ./ d')  # no need to update AA, only B is needed.

        ng = norm(g)
        if verbose
            println("Iter $iter crit = $o1 ng = $ng")
        end
        count_iter += 1
        if ng < tol * o1
            break
        end
    end

    # now, we can update AA
    for q = 1:Q
        AA[:, :, q] = L .* AA[:, :, q] ./ L'
    end

    return L, AA, count_iter
end
