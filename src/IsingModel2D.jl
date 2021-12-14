module IsingModel2D
    using Statistics, LinearAlgebra, Random
    import SpecialFunctions

    include("onsager.jl")
    include("ising.jl")

    include("metropolis.jl")
    include("metropolis_f.jl")
    include("wolff.jl")
    include("hybrid.jl")
end # module
