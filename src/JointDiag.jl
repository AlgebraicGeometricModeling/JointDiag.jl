module JointDiag

using LinearAlgebra

import Random

"""
    AbstractSolver

Solver for joint diagonalization of matrices.
"""
abstract type AbstractSolver end

function joint_diag(
    matrices::AbstractVector{<:AbstractMatrix{T}},
    solver::AbstractSolver,
) where {T}
    λ = rand(float(T), length(matrices))
    λ ./= sum(λ)
    return joint_diag(matrices, λ, solver)
end

function joint_diag(
    matrices::AbstractVector{<:AbstractMatrix},
    λ::AbstractVector,
    solver::AbstractSolver,
)
    @assert length(matrices) == length(λ)
    return joint_diag(matrices, sum(λ .* matrices), solver)
end

function solution_matrix(
    matrices::AbstractVector{<:AbstractMatrix{T}},
    left_factor,
    right_factor,
) where {T}
    X = zeros(eltype(left_factor), length(matrices), size(right_factor, 2))
    for j in axes(X, 2)
        left = left_factor[j, :]
        right = right_factor[:, j]
        for i in axes(X, 1)
            X[i, j] = LinearAlgebra.transpose(left) * matrices[i] * right
        end
    end
    return X
end

include("joint_diag_rand.jl")
include("joint_diag_newton.jl")

include("cluster.jl")
include("schur.jl")

include("joint_reduce.jl")
end
