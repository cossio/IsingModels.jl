using Test: @test, @testset, @inferred
using Random: bitrand
import IsingModels as Ising

@testset "wolff Padd" begin
    for β in 0:0.1:1
        @test Ising.wolff_padd(β) ≈ 1 - exp(-2β)
    end
end

@testset "wolff" begin
    L = 50; T = 100; Δ = 12
    σ = bitrand(L, L)
    σ_t, M, E = @inferred Ising.wolff!(σ, Ising.βc; steps=T, save_interval=Δ)
    @test length(M) == length(E) == T
    @test size(σ_t) == (size(σ)..., length(1:Δ:T))
    @test M[end] == Ising.magnetization(σ)
    @test E[end] ≈ Ising.energy(σ)
    @test M[1:Δ:end] == Ising.magnetization.(eachslice(σ_t; dims=3))
    @test E[1:Δ:end] == Ising.energy.(eachslice(σ_t; dims=3))
end

@testset "wolff energy at every step" begin
    # dimensions of size 1 and 2 have coinciding neighbors; also 1-D and 3-D lattices
    sizes = ((7,), (2, 5), (1, 6), (8, 8), (2, 2, 2), (3, 2, 4), (5, 5, 5))
    for sz in sizes, β in (0, 0.2, Ising.βc, 1, Inf)
        σ = bitrand(sz)
        σ_t, M, E = Ising.wolff!(σ, β; steps=300, save_interval=1)
        @test M == Ising.magnetization.(eachslice(σ_t; dims=length(sz) + 1))
        @test E == Ising.energy.(eachslice(σ_t; dims=length(sz) + 1))
        @test E[end] == Ising.energy(σ)
    end
end

@testset "wolff_step! cluster buffer" begin
    σ = bitrand(20, 20)
    σ0 = copy(σ)
    M = [Ising.magnetization(σ), 0]
    E = [Ising.energy(σ), 0]
    cluster = falses(size(σ))
    cluster_size = Ising.wolff_step!(σ; t=2, M, E, Padd=Ising.wolff_padd(Ising.βc), cluster)
    @test !any(cluster) # left clean for the next step
    @test count(σ .!= σ0) == cluster_size
    @test M[2] == Ising.magnetization(σ)
    @test E[2] == Ising.energy(σ)
end

@testset "wolff cluster" begin
    σ = bitrand(50, 50)
    @inferred Ising.metropolis!(σ, Ising.βc; steps=10^6)
    i_seed = rand(CartesianIndices(σ))
    cluster = @inferred Ising.wolff_cluster(σ, i_seed, Ising.wolff_padd(Ising.βc))

    # seed is always flipped
    @test cluster[i_seed]

    # all cluster spins have same sign as the center
    @test cluster .* σ == σ[i_seed] * cluster

    # wolff_flip! marks the sites it returns, flips them, and returns the energy change
    σ_flip = copy(σ)
    cluster = falses(size(σ))
    sites, ΔE = Ising.wolff_flip!(cluster, σ_flip, i_seed, Ising.wolff_padd(Ising.βc))
    @test sort(sites) == findall(cluster)
    @test (σ_flip .!= σ) == cluster
    @test ΔE == Ising.energy(σ_flip) - Ising.energy(σ)

    # zero temperature cluster
    σ = bitrand(50, 50)
    clset = Set([CartesianIndex(1,1)])
    queue = [first(clset)]
    while !isempty(queue)
        i = pop!(queue)
        for x in Ising.neighbors(CartesianIndices(σ), i)
            if σ[x] == σ[i] && x ∉ clset
                push!(clset, x)
                push!(queue, x)
            end
        end
    end
    cluster = CartesianIndices(σ) .∈ Ref(clset)
    #cluster = BitMatrix([(i,j) ∈ cluster for i in 1:50, j in 1:50])
    @test cluster[CartesianIndex(1,1)]
    @test Ising.wolff_cluster(σ, CartesianIndex(1,1)) == cluster

    cluster = Ising.wolff_cluster(σ, CartesianIndex(1,1))
    for i in CartesianIndices(σ)
        if cluster[i]
            for x in Ising.neighbors(CartesianIndices(σ), i)
                if !cluster[x]
                    # at zero temperature, spins not in cluster have different sign
                    @test σ[x] ≠ σ[CartesianIndex(1,1)]
                end
            end
        end
    end
end

@testset "flip_cluster" begin
    σ = falses(64, 64)
    c = bitrand(64, 64)
    σ0 = copy(σ)
    Ising.flip_cluster!(σ, c)
    @test (σ0 .!= σ) == c
end
