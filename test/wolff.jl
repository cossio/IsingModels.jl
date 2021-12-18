include("init.jl")

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

@testset "wolff cluster" begin
    σ = bitrand(50, 50)
    @inferred Ising.metropolis!(σ, Ising.βc; steps=10^6)
    i_seed = rand(CartesianIndices(σ))
    cluster = @inferred Ising.wolff_cluster(σ, i_seed, Ising.wolff_padd(Ising.βc))

    # seed is always flipped
    @test cluster[i_seed]

    # all cluster spins have same sign as the center
    @test cluster .* σ == σ[i_seed] * cluster

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
