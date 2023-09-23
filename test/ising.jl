using Test: @test, @testset, @inferred
using Random: bitrand
import IsingModels as Ising

@test Ising.spin(true) == 1
@test Ising.spin(false) == -1
@test Ising.binary(-2) == false
@test Ising.binary(0)  == false
@test Ising.binary(2) == true

for σ in (true, false)
    @test Ising.binary(Ising.spin(σ)) == σ
end

for s in (-1, 1)
    @test Ising.spin(Ising.binary(s)) == s
end

@testset "neighbors" begin
    σ = trues(3, 3)
    for i in CartesianIndices(σ)
        @test Ising.sum_neighbors(σ, i) == 4
    end

    σ = falses(3, 3)
    for i in CartesianIndices(σ)
        @test Ising.sum_neighbors(σ, i) == -4
    end

    σ = bitrand(3, 3)
    Ising.sum_neighbors(σ, CartesianIndex(2,2)) == (
        Ising.spin(σ[1,2]) + Ising.spin(σ[3,2]) +
        Ising.spin(σ[2,1]) + Ising.spin(σ[2,3])
    )

    for D in (2, 3)
        σ = bitrand(rand(32:64, D)...)
        for i in CartesianIndices(σ)
            @test Ising.sum_neighbors(σ, i) ∈ -2D:2:2D
        end
        M = Ising.magnetization(σ)
        @test sum(Ising.sum_neighbors(σ, i) for i in CartesianIndices(σ)) == 2D * M
    end
end

@testset "Ising 1-D" begin
    L = 1000
    σ = bitrand(1000)
    for i in 1:L
        S = Ising.spin(σ[mod1(i - 1, L)]) + Ising.spin(σ[mod1(i + 1, L)])
        @test Ising.sum_neighbors(σ, CartesianIndex(i)) == S
    end
end

@testset "magnetization" begin
    L, K = 10, 15
    σ = bitrand(L, K)
    @test Ising.magnetization(σ) == sum(Ising.spin, σ)
end

@testset "energy" begin
    L, K = 10, 15
    σ = bitrand(L, K)
    E = -sum(Ising.spin(σ[i]) * Ising.sum_neighbors(σ,i) for i in CartesianIndices(σ))
    @test 2Ising.energy(σ) == E
    h = randn()
    @test Ising.energy(σ, h) ≈ Ising.energy(σ) - h * Ising.magnetization(σ)
end
