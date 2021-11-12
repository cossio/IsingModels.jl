include("init.jl")

@testset "Ising grid lattice" begin
    @test Ising.βc ≈ 0.44068679350977151261630466248989615451408016413081770537664780432668
    @test iszero(Ising.onsager_magnetization(Ising.βc))
    @test iszero(Ising.onsager_magnetization(0))
    @test iszero(Ising.onsager_magnetization(0.1))
    @test Ising.onsager_magnetization(100) ≈ 1

    spins = rand((Int8(-1), Int8(1)), 50, 50)
    @test neighbor_sum(spins, 5, 5) == spins[5,4] + spins[5,6] + spins[4,5] + spins[6,5]
end
