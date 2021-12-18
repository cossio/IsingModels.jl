include("init.jl")

@testset "metropolis" begin
    L = 50; T = 100; Δ = 12; β = 1.0
    σ = bitrand(L, L)
    σ_t, M, E = @inferred Ising.metropolis!(σ, β; steps=T, save_interval=Δ)
    @test length(M) == length(E) == T
    @test size(σ_t) == (size(σ)..., length(1:Δ:T))
    @test M[end] == @inferred Ising.magnetization(σ)
    @test E[end] ≈ @inferred Ising.energy(σ)
    @test M[1:Δ:end] == Ising.magnetization.(eachslice(σ_t; dims=3))
    @test E[1:Δ:end] == Ising.energy.(eachslice(σ_t; dims=3))
end

@testset "metropolis + f(M)" begin
    L = 50; T = 100; Δ = 12; β = 1.0
    f = sin
    σ = bitrand(L, L)
    σ_t, M, E = @inferred Ising.metropolis!(σ, β; steps=T, save_interval=Δ, f=f)
    @test length(M) == length(E) == T
    @test size(σ_t) == (size(σ)..., length(1:Δ:T))
    @test M[end] ≈ @inferred Ising.magnetization(σ)
    @test E[end] ≈ @inferred Ising.energy(σ; f=f)
    @test M[1:Δ:end] ≈ Ising.magnetization.(eachslice(σ_t; dims=3))
    @test E[1:Δ:end] ≈ Ising.energy.(eachslice(σ_t; dims=3); f=f)
end
