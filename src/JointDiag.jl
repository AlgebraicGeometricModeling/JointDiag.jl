module JointDiag

using LinearAlgebra
import Random

"""
    AbstractSolver

Solver for joint diagonalization of matrices.
"""
abstract type AbstractSolver end

include("diagonalisation.jl")
include("cluster.jl")
include("schur.jl")

end
