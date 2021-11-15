include("init.jl")

@testset "metropolis" begin
    spins = Ising.random_configuration(50)
    spins_t, M, E = @inferred Ising.metropolis!(spins, 1.0, 100)
    @test M[end] == sum(spins)
    @test E[end] ≈ @inferred Ising.energy(spins)
end
