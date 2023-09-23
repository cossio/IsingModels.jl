using Test: @testset, @test, @inferred
using Random: bitrand
import IsingModels as Ising

@testset "Curie Weiss model" begin
    N = 128
    σ = bitrand(N)
    @test length(σ) == N

    M = sum(Ising.spin, σ)
    E = -sum(Ising.spin(σ[i]) * Ising.spin(σ[j]) for i = 1:N for j = 1:(i - 1)) / N

    @test M ≈ @inferred Ising.magnetization(σ)
    @test E ≈ @inferred Ising.mf_energy(σ)

    @test E - π * M ≈ @inferred Ising.mf_energy(σ, π)
    @test E + sin(M) ≈ @inferred Ising.mf_energy(σ; f=sin)
    @test E - π * M + sin(M) ≈ @inferred Ising.mf_energy(σ, π; f=sin)

    @test all(Ising.mf_random_magnetized_configuration(N, N))
    @test !any(Ising.mf_random_magnetized_configuration(N, -N))
    M = 20
    σ = Ising.mf_random_magnetized_configuration(N, M)
    @test Ising.magnetization(σ) == M
end
