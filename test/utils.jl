function testelements(X, Y; atol = Base.rtoldefault(Float64), kwargs...)
    @test size(X) == size(Y)
    for y in eachcol(Y)
        @test any(x -> isapprox(x, y; atol = atol, kwargs...), eachcol(X))
    end
end

function testelements(Y, args...; kwargs...)
    X, _, solver = JointDiag.joint_diag(args...)
    Xc = JointDiag.cluster(X, solver, Base.rtoldefault(Float64))
    testelements(Xc, Y; kwargs...)
end

function testelementstypes(X, Y; kwargs...)
    testelements(X, Y; kwargs...)
    for T in [Rational{Int}, Float64]
        if X isa FixedVariablesSet
            U = T
        else
            U = float(T)
        end
        V = similar(X, T)
        testelements(V, Y; kwargs...)
        @test eltype(V) == Vector{U}
        @test collect(V) isa Vector{Vector{U}}
    end
end

# We use a fixed RNG in the tests to decrease nondeterminism. There is still nondeterminism in LAPACK though
schur_solver = JointDiag.SchurJointDiag()
