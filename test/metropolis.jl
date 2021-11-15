include("init.jl")

@testset "metropolis" begin
    spins = Ising.random_configuration(50)
    spins_t, m, E = @inferred Ising.metropolis!(spins, 1.0, 100)
    @test round(Int, m[end] * length(spins)) == sum(spins)
    @test E[end] ≈ @inferred Ising.energy(spins) / length(spins)
end
