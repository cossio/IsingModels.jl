using Test: @test, @testset, @inferred
using Random: bitrand
using Zygote
import IsingModels as Ising

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
    for β = 0.1:0.1:1
        @test -β^2 * Ising.onsager_internal_energy'(β) ≈ Ising.onsager_heat_capacity(β)
    end
end

@testset "Kramers-Wannier duality of internal energy and heat capacity" begin
    for β in 0.1:0.1:1
        β_dual = Ising.kramers_wannier(β)
        E_here = Ising.onsager_internal_energy(β)
        E_dual = Ising.onsager_internal_energy(β_dual)
        @test E_dual ≈ -2cosh(2β) - sinh(2β) * E_here

        C_here = Ising.onsager_heat_capacity(β)
        C_dual = Ising.onsager_heat_capacity(β_dual)

        x_dual = C_dual / β_dual^2
        x_here = C_here / β^2
        @test x_dual ≈ (x_here - 4) * sinh(2β)^2 - sinh(4β) * E_here
    end
end
