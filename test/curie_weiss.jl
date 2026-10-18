using Test: @testset, @test, @inferred
using Random: bitrand
using Statistics: mean
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

@testset "Curie Weiss sampling" begin
    N = 10; β = 0.7; h = 0.2; B = 10^5
    f(M) = M^2 / 100

    # exact distribution of M, by enumeration of all configurations
    Ms = -N:2:N
    P = zeros(length(Ms))
    for k in 0:(2^N - 1)
        σ = BitVector(digits(k; base=2, pad=N) .== 1)
        P[(Ising.magnetization(σ) + N) ÷ 2 + 1] += exp(-β * Ising.mf_energy(σ, h; f))
    end
    P /= sum(P)

    samples = Ising.mf_random_magnetizations(N, β, h; f, B)
    @test length(samples) == B
    @test all(∈(Ms), samples)
    @test [mean(==(M), samples) for M in Ms] ≈ P atol=0.01

    # default f = nothing and h = 0, which is symmetric under M → -M
    samples = Ising.mf_random_magnetizations(N, β; B)
    @test all(∈(Ms), samples)
    @test abs(mean(samples)) < 0.1

    σ = Ising.mf_random_configurations(N, β, h; f, B=7)
    @test size(σ) == (N, 7)
    σ = Ising.mf_random_configurations(N, β)
    @test size(σ) == (N, 1)
end
