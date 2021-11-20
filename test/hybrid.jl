include("init.jl")

@testset "hybrid" begin
    L = 50; T = 10000; Δsave = 12; Δlocal = L; β = Ising.βc
    spins = Ising.random_configuration(L)
    spins_t, M, E = @inferred Ising.hybrid!(spins, β, T; save_interval=Δsave, local_steps=Δlocal)
    @test length(M) == length(E) == T
    @test size(spins_t) == (size(spins)..., length(1:Δsave:T))
    @test M[end] == sum(spins)
    @test E[end] ≈ Ising.energy(spins)
    @test M[1:Δsave:end] == dropdims(sum(spins_t; dims=(1,2)); dims=(1,2))
    @test E[1:Δsave:end] == Ising.energy.(eachslice(spins_t; dims=3))
end

@testset "dynamic_hybrid" begin
    L = 50; T = 10000; Δsave = 12; β = Ising.βc
    spins = Ising.random_configuration(L)
    spins_t, M, E = @inferred Ising.dynamic_hybrid!(spins, β, T; save_interval=Δsave)
    @test length(M) == length(E) == T
    @test size(spins_t) == (size(spins)..., length(1:Δsave:T))
    @test M[end] == sum(spins)
    @test E[end] ≈ Ising.energy(spins)
    @test M[1:Δsave:end] == dropdims(sum(spins_t; dims=(1,2)); dims=(1,2))
    @test E[1:Δsave:end] == Ising.energy.(eachslice(spins_t; dims=3))
end
