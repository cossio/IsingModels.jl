module IsingModels
    using Statistics, LinearAlgebra, Random
    import SpecialFunctions

    include("onsager.jl")
    include("ising.jl")

    include("metropolis.jl")
    include("metropolis_f.jl")
    include("wolff.jl")
    include("hybrid.jl")

    include("curie_weiss.jl")
end # module
