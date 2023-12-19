module QPN

using CDDLib
using Combinatorics
using DataStructures
using IterTools
using LinearAlgebra
using OSQP
using PATHSolver
using Polyhedra
using Printf
using ProgressMeter
using Random
using SparseArrays
using Symbolics
using GLMakie

using Infiltrator

#PATHSolver.c_api_License_SetString("2830898829&Courtesy&&&USR&45321&5_1_2021&1000&PATH&GEN&31_12_2025&0_0_0&6000&0_0")

include("macros.jl")
include("sets.jl")
include("intersection.jl")
include("programs.jl")
include("algorithm.jl")
include("qp_processing.jl")
include("requests.jl")
include("avi2.jl")
include("avi_solutions2.jl")
include("printing.jl")

include("../examples/bilevel_escape.jl")
include("../examples/trilevel_escape.jl")
include("../examples/repeated_variable_control.jl")
include("../examples/simple_bilevel.jl")
include("../examples/simple_network.jl")
include("../examples/robust_constrained.jl")
include("../examples/control_avoid.jl")
include("../examples/robust_avoid.jl")
include("../examples/rock_paper_scissors.jl")

export Poly, QP, Constraint, Quadratic, QPNet, solve, setup

end # module QPN
