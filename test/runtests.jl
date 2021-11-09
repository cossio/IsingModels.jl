using Test, SafeTestsets, Aqua, SquareIsingModel
include("init.jl")

@safetestset "util" begin include("util.jl") end
@safetestset "ancestors" begin include("ancestors.jl") end
@safetestset "energies" begin include("energies.jl") end
@safetestset "model" begin include("model.jl") end
@safetestset "data" begin include("data.jl") end
@safetestset "rbm" begin include("rbm.jl") end
@safetestset "tree" begin include("tree.jl") end

#Aqua.test_all(PhageTree)
