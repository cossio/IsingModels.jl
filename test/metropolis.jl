include("init.jl")

@testset "metropolis" begin
    spins = rand((Int8(-1), Int8(1)), 50, 50)
    @test neighbor_sum(spins, 5, 5) == spins[5,4] + spins[5,6] + spins[4,5] + spins[6,5]
    spins_t, m, E = Ising.metropolis!(spins, 1.0, 100)
    @test m[end] ≈ mean(spins)
    @test E[end] ≈ Ising.energy(spins) / length(spins)
end
