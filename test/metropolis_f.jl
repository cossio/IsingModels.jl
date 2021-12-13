include("init.jl")

@testset "metropolis_f" begin
    L = 50; T = 100; Δ = 12; β = 1.0
    f = sin
    spins = Ising.random_configuration(L)
    spins_t, M, E = @inferred Ising.metropolis_f!(spins, β; steps=T, save_interval=Δ, f=f)
    @test length(M) == length(E) == T
    @test size(spins_t) == (size(spins)..., length(1:Δ:T))
    @test M[end] == Ising.magnetization(spins)
    @test E[end] ≈ @inferred Ising.energy(spins) + f(M[end])
    @test M[1:Δ:end] == Ising.magnetization.(eachslice(spins_t; dims=3))
    @test E[1:Δ:end] ≈ Ising.energy.(eachslice(spins_t; dims=3)) + f.(M[1:Δ:end])
end
