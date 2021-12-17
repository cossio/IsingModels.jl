using Test, SafeTestsets
import IsingModels as Ising

@safetestset "onsager" begin include("onsager.jl") end
@safetestset "ising" begin include("ising.jl") end
#@safetestset "curie_weiss" begin include("curie_weiss.jl") end
@safetestset "metropolis" begin include("metropolis.jl") end
@safetestset "metropolis_f" begin include("metropolis_f.jl") end
@safetestset "wolff" begin include("wolff.jl") end
@safetestset "hybrid" begin include("hybrid.jl") end
@safetestset "bit" begin include("bit.jl") end
