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

include("joint_diag_rand.jl")
include("joint_diag_newton.jl")

include("cluster.jl")
include("schur.jl")

include("joint_reduce.jl")
end
