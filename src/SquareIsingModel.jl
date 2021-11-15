module SquareIsingModel
    using Statistics, LinearAlgebra

    include("onsager.jl")
    include("ising.jl")

    include("metropolis.jl")
    include("wolff.jl")
    include("hybrid.jl")
end # module
