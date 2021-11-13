include("init.jl")

@testset "metropolis" begin
    spins = Ising.random_configuration(50, 50)
    spins_t, m, E = Ising.metropolis!(spins, 1.0, 100)
    @test m[end] ≈ mean(spins)
    @test E[end] ≈ Ising.energy(spins) / length(spins)
end
