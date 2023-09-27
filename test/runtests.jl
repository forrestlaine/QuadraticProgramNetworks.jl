using Test
using QPN

# Example tests
include("repeated_variable_control.jl")
include("trilevel_escape.jl")
include("bilevel_escape.jl")
include("simple_bilevel.jl")

# Low-level tests
include("convex_hull.jl")
