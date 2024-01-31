module JointDiag

using LinearAlgebra
import Random

"""
    AbstractSolver

Solver for joint diagonalization of matrices.
"""
abstract type AbstractSolver end

include("rand_joint_diag.jl")
include("newton_joint_diag.jl")
include("cluster.jl")
include("schur.jl")

end
