include("init.jl")

@testset "Ising grid lattice" begin
    @test Ising.βc ≈ 0.44068679350977151261630466248989615451408016413081770537664780432668
    @test iszero(Ising.onsager_magnetization(Ising.βc))
    @test iszero(Ising.onsager_magnetization(0))
    @test iszero(Ising.onsager_magnetization(0.1))
    @test Ising.onsager_magnetization(100) ≈ 1

    s = @inferred Ising.random_configuration(50)
    @test Ising.neighbor_sum(s, 5, 5) == s[5,4] + s[5,6] + s[4,5] + s[6,5]

    @test size(Ising.random_configuration(50)) == (50, 50)
    @test size(Ising.random_configuration(20, 40)) == (20, 40)
end
