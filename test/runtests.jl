using Test, SafeTestsets, SquareIsingModel

@safetestset "onsager" begin include("onsager.jl") end
@safetestset "ising" begin include("ising.jl") end
@safetestset "metropolis" begin include("metropolis.jl") end
@safetestset "wolff" begin include("wolff.jl") end
