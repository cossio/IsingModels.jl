using Test: @test, @testset, @inferred
using Random: bitrand
using Statistics: mean
import IsingModels as Ising

@testset "hybrid" begin
    L = 50; T = 10000; Δsave = 12; Δlocal = L
    σ = bitrand(L, L)
    σ_t, M, E = @inferred Ising.hybrid!(σ, Ising.βc; steps=T, save_interval=Δsave, local_steps=Δlocal)
    @test length(M) == length(E) == T
    @test size(σ_t) == (size(σ)..., length(1:Δsave:T))
    @test M[end] == Ising.magnetization(σ)
    @test E[end] ≈ Ising.energy(σ)
    @test M[1:Δsave:end] == Ising.magnetization.(eachslice(σ_t; dims=3))
    @test E[1:Δsave:end] == Ising.energy.(eachslice(σ_t; dims=3))
end

@testset "dynamic_hybrid" begin
    L = 50; T = 10000; Δsave = 12
    σ = bitrand(L, L)
    σ_t, M, E = @inferred Ising.dynamic_hybrid!(σ, Ising.βc; steps=T, save_interval=Δsave)
    @test length(M) == length(E) == T
    @test size(σ_t) == (size(σ)..., length(1:Δsave:T))
    @test M[end] == Ising.magnetization(σ)
    @test E[end] ≈ Ising.energy(σ)
    @test M[1:Δsave:end] == Ising.magnetization.(eachslice(σ_t; dims=3))
    @test E[1:Δsave:end] == Ising.energy.(eachslice(σ_t; dims=3))
end

@testset "hybrid, D = $D" for D in (1, 3)
    T = 1000; Δsave = 12
    σ = bitrand(ntuple(_ -> 6, D)...)
    for sampler! in (Ising.hybrid!, Ising.dynamic_hybrid!)
        σ_t, M, E = sampler!(σ, Ising.βc; steps=T, save_interval=Δsave)
        @test size(σ_t) == (size(σ)..., length(1:Δsave:T))
        @test M[end] == Ising.magnetization(σ)
        @test E[end] == Ising.energy(σ)
        @test M[1:Δsave:end] == Ising.magnetization.(eachslice(σ_t; dims=D + 1))
        @test E[1:Δsave:end] == Ising.energy.(eachslice(σ_t; dims=D + 1))
    end
end

@testset "hybrid samplers vs. exact enumeration" begin
    # 3x3 lattice, small enough to enumerate all configurations
    L = 3; β = 0.4
    confs = [BitMatrix(reshape(digits(k; base=2, pad=L^2) .== 1, L, L)) for k in 0:(2^(L^2) - 1)]
    Es = Ising.energy.(confs)
    w = exp.(-β * Es)
    E_exact = sum(w .* Es) / sum(w)
    for sampler! in (Ising.hybrid!, Ising.dynamic_hybrid!)
        _, _, E = sampler!(bitrand(L, L), β; steps=10^6)
        @test mean(E) ≈ E_exact atol=0.3
    end
end
