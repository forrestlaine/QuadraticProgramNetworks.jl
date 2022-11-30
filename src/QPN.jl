module QPN

using Combinatorics
using CDDLib
using DataStructures
using LinearAlgebra
using IterTools
using OSQP
using PATHSolver
using Printf
using ProgressMeter
using Polyhedra
using Random
using SparseArrays

using Infiltrator

PATHSolver.c_api_License_SetString("2830898829&Courtesy&&&USR&45321&5_1_2021&1000&PATH&GEN&31_12_2025&0_0_0&6000&0_0")

include("macros.jl")
include("sets.jl")
include("intersection.jl")
include("programs.jl")
include("algorithm.jl")
include("avi.jl")
include("avi_solutions.jl")
include("printing.jl")

export Poly, QP, Quadratic, QPNet, solve

end # module QPN
