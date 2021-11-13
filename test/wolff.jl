include("init.jl")

@testset "wolff" begin
    spins = Ising.random_configuration(50, 50)
    spins_t, m, E = Ising.wolff!(spins, 1.0, 10)
    @test m[end] ≈ mean(spins)
    @test E[end] ≈ Ising.energy(spins) / length(spins)
end

@testset "wolff cluster" begin
    spins = Ising.random_configuration(50, 50)
    β = 0.5
    Ising.metropolis!(spins, β, 10^6)
    cluster = Ising.wolff_cluster(spins, 10, 10, 1 - exp(-2β))
    @test cluster .* spins == spins[10,10] * cluster
end
