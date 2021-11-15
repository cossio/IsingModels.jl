include("init.jl")

@testset "metropolis acceptance probabilities" begin
    β = Ising.βc
    _Exp2β = Ising.metropolis_acceptance_probabilities(β)
    @test length(_Exp2β) == 5
end

@testset "metropolis" begin
    L = 50; T = 100; Δ = 12; β = 1.0
    spins = Ising.random_configuration(L)
    spins_t, M, E = @inferred Ising.metropolis!(spins, β, T; save_interval=Δ)
    @test length(M) == length(E) == T
    @test size(spins_t) == (size(spins)..., length(1:Δ:T))
    @test M[end] == sum(spins)
    @test E[end] ≈ @inferred Ising.energy(spins)
    @test M[1:Δ:end] == dropdims(sum(spins_t; dims=(1,2)); dims=(1,2))
    @test E[1:Δ:end] == Ising.energy.(eachslice(spins_t; dims=3))
end
