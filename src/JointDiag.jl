module JointDiag

using LinearAlgebra

import Random

"""
    AbstractSolver

Solver for joint diagonalization of matrices.
"""
abstract type AbstractSolver end

include("joint_diag_rand.jl")
include("joint_diag_newton.jl")
include("joint_diag_jacobi.jl")

include("cluster.jl")
include("schur.jl")

include("joint_reduce.jl")
end
