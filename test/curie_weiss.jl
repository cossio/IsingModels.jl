include("init.jl")

@testset "Curie Weiss model" begin
    N = 100; T = 100; Δ = 12; β = 1.0
    model = Ising.CurieWeiss(N)
    @test length(model.spins) == N

    M = sum(model.spins)
    @test M ≈ @inferred Ising.magnetization(model)

    E = -sum(model.spins[i] * model.spins[j] for i = 1:N for j = 1:(i - 1)) / N
    @test E - π * M ≈ @inferred Ising.energy(model, π)
    @test E + sin(M) ≈ @inferred Ising.energy(model, f=sin)

    spins_t, M, E = @inferred Ising.metropolis!(model, β; steps=T, save_interval=Δ)
    @test length(M) == length(E) == T
    @test size(spins_t) == (N, length(1:Δ:T))
    @test M[end] ≈ @inferred Ising.magnetization(model)
    @test E[end] ≈ @inferred Ising.energy(model)
    @test M[1:Δ:end] ≈ Ising.magnetization.(Ising.CurieWeiss.(eachcol(spins_t)))
    @test E[1:Δ:end] ≈ Ising.energy.(Ising.CurieWeiss.(eachcol(spins_t)))
end
