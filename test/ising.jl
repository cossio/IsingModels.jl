include("init.jl")

@testset "Ising grid lattice" begin
    s = @inferred Ising.random_configuration(50)
    @test Ising.neighbor_sum(s, 5, 5) == s[5,4] + s[5,6] + s[4,5] + s[6,5]
    @test size(Ising.random_configuration(50)) == (50, 50)
    @test size(Ising.random_configuration(20, 40)) == (20, 40)
end

@testset "neighbor sum" begin
    spins = Int8[
        +1  +1  +1;
        +1  +1  +1;
        +1  +1  +1;
    ]
    @test Ising.neighbor_sum(spins, 2, 2) == 4

    spins = Int8[
        +1  -1  +1;
        +1  +1  +1;
        +1  +1  +1;
    ]
    @test Ising.neighbor_sum(spins, 2, 2) == 2

    spins = Int8[
        +1  -1  +1;
        +1  +1  -1;
        +1  +1  +1;
    ]
    @test Ising.neighbor_sum(spins, 2, 2) == 0

    spins = Int8[
        +1  -1  +1;
        +1  +1  -1;
        +1  -1  +1;
    ]
    @test Ising.neighbor_sum(spins, 2, 2) == -2

    spins = Int8[
        +1  -1  +1;
        -1  +1  -1;
        +1  -1  +1;
    ]
    @test Ising.neighbor_sum(spins, 2, 2) == -4

    L = 100
    spins = Ising.random_configuration(L)
    for i = 1:L, j = 1:L
        @test Ising.neighbor_sum(spins, i, j) ∈ -4:2:4
    end
end

@testset "distance tensor" begin
    L = 10
    d = Ising.distance_tensor(L)
    for i1 = 1:L, j1 = 1:L, i2 = 1:L, j2 = 1:L
        di = min(abs(i1 - i2), abs(i1 - i2 - L), abs(i1 - i2 + L))
        dj = min(abs(j1 - j2), abs(j1 - j2 - L), abs(j1 - j2 + L))
        @test d[i1,j1,i2,j2] ≈ sqrt(di^2 + dj^2)
    end
end
