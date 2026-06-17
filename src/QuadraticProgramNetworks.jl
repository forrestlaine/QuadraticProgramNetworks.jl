module QuadraticProgramNetworks

using CDDLib
using Combinatorics
using LinearAlgebra
using OSQP
using PATHSolver
using Polyhedra
using Printf
using ProgressMeter
using Random
using SparseArrays
using Symbolics

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

export Poly,
       QP,
       Constraint,
       Quadratic,
       QPNet,
       solve,
       setup,
       variable,
       variables,
       add_constraint!,
       add_qp!,
       add_edges!,
       assign_constraint_groups!,
       set_options!,
       display_solution

end # module QPN
