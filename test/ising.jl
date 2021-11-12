include("init.jl")

@testset "Ising grid lattice" begin
    @test Ising.βc ≈ 0.44068679350977151261630466248989615451408016413081770537664780432668
    @test iszero(Ising.onsager_magnetization(Ising.βc))
    @test iszero(Ising.onsager_magnetization(0))
    @test iszero(Ising.onsager_magnetization(0.1))
    @test Ising.onsager_magnetization(100) ≈ 1
end
