include("init.jl")

@test Ising.βc ≈ 0.44068679350977151261630466248989615451408016413081770537664780432668

@testset "Onsager magnetization" begin
    @test iszero(Ising.onsager_magnetization(Ising.βc))
    @test iszero(Ising.onsager_magnetization(0))
    @test iszero(Ising.onsager_magnetization(0.1))
    @test Ising.onsager_magnetization(100) ≈ 1
end

@testset "Kramers-Wannier duality" begin
    @test Ising.kramers_wannier(Ising.βc) ≈ Ising.βc
    @test isinf(Ising.kramers_wannier(0))
    @test iszero(Ising.kramers_wannier(Inf))
    for β in 0.1:0.1:1
        β_ = Ising.kramers_wannier(β)
        @test sinh(2β_) * sinh(2β) ≈ 1
    end
end

@testset "Onsager heat capacity and energy thermodynamic relations" begin
    u(T) = Ising.onsager_internal_energy(1/T)
    for T in 0.5:0.5:5
        @test u'(T) ≈ Ising.onsager_heat_capacity(1/T) rtol=0.0001
    end
end

@testset "Kramers-Wannier duality of heat capacity and internal energy" begin
    for β in 0.1:0.1:1
        Ehere = Ising.onsager_internal_energy(β)
        Edual = Ising.onsager_internal_energy(Ising.kramers_wannier(β))
        @test Edual ≈ -2cosh(2β) - sinh(2β) * Ehere
    end
end
