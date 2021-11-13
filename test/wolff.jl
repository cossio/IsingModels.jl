include("init.jl")

@testset "wolff Padd" begin
    for β in 0:0.1:1
        @test Ising.wolff_padd(β) ≈ 1 - exp(-2β)
    end
end

@testset "wolff" begin
    spins = Ising.random_configuration(50)
    spins_t, m, E = @inferred Ising.wolff!(spins, 1.0, 10)
    @test m[end] ≈ mean(spins)
    @test E[end] ≈ Ising.energy(spins) / length(spins)
end

@testset "wolff cluster" begin
    spins = Ising.random_configuration(50)
    β = 0.5
    @inferred Ising.metropolis!(spins, β, 10^6)
    cluster = @inferred Ising.wolff_cluster(spins, 1, 1, Ising.wolff_padd(β))
    @test cluster[1, 1]
    # all cluster spins have same sign as the center
    @test cluster .* spins == spins[1,1] * cluster

    # zero temperature cluster
    spins = Ising.random_configuration(50)
    cluster = Set([(1,1)])
    queue = [(1,1)]
    while !isempty(queue)
        (i,j) = pop!(queue)
        for (x,y) in Ising.neighbors(i,j, size(spins)...)
            if spins[x,y] == spins[i,j] && (x,y) ∉ cluster
                push!(cluster, (x,y))
                push!(queue, (x,y))
            end
        end
    end
    cluster = BitMatrix([(i,j) ∈ cluster for i in 1:50, j in 1:50])
    @test Ising.wolff_cluster(spins, 1, 1) == cluster
    @test cluster[1,1]

    cluster = Ising.wolff_cluster(spins, 1, 1)
    for i in 1:50, j in 1:50
        if cluster[i,j]
            for (x, y) in Ising.neighbors(i, j, 50, 50)
                if !cluster[x,y]
                    # at zero temperature, spins not in cluster have different sign
                    @test spins[x,y] ≠ spins[1, 1]
                end
            end
        end
    end
end
