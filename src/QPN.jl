module QPN

using CDDLib
using Combinatorics
using DataStructures
using IterTools
using LinearAlgebra
using OSQP
using PATHSolver
using LaTeXStrings
using Polyhedra
using Printf
using ProgressMeter
using Random
using SparseArrays
using Symbolics
using CairoMakie

using Infiltrator

include("macros.jl")
include("sets.jl")
include("intersection.jl")
include("programs.jl")
include("algorithm.jl")
include("requests.jl")
include("qp_processing.jl")
include("avi.jl")
include("avi_solutions.jl")
include("printing.jl")

include("../examples/bilevel_escape.jl")
include("../examples/trilevel_escape.jl")
include("../examples/repeated_variable_control.jl")
include("../examples/simple_bilevel.jl")
include("../examples/simple_network.jl")
include("../examples/robust_constrained.jl")
include("../examples/control_avoid.jl")
include("../examples/robust_avoid.jl")
include("../examples/robust_avoid_simple.jl")
include("../examples/four_player_matrix_game.jl")
include("../examples/rock_paper_scissors.jl")

export Poly, QP, Constraint, Quadratic, QPNet, solve, setup

end # module QPN
