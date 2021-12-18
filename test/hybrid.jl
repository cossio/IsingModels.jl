include("init.jl")

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
