include("init.jl")

@testset "wolff" begin
    spins = rand((Int8(-1), Int8(1)), 50, 50)
    spins_t, m, E = Ising.wolff!(spins, 1.0, 10)
    @test m[end] ≈ mean(spins)
    @test E[end] ≈ Ising.energy(spins) / length(spins)
end
