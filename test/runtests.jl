using Test
using QuadraticProgramNetworks

@testset "package boundary" begin
    @test !isdefined(QuadraticProgramNetworks, :visualize_simple_bilevel)
    @test !isdefined(QuadraticProgramNetworks, :viz_solution_robust_avoid_simple)
end

include("simple_bilevel.jl")
