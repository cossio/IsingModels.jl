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

@testset "Onsager internal energy at criticality" begin
    @test Ising.onsager_internal_energy(Ising.βc) ≈ -√2
    @test Ising.onsager_internal_energy(prevfloat(Ising.βc)) ≈ -√2
    @test Ising.onsager_internal_energy(nextfloat(Ising.βc)) ≈ -√2
    @test Ising.onsager_internal_energy(Ising.βc + 1e-9) ≈ -√2 rtol=1e-6
    @test Ising.onsager_internal_energy(Ising.βc - 1e-9) ≈ -√2 rtol=1e-6
    for β in Ising.βc .+ (-1e-7:1e-9:1e-7)
        @test isfinite(Ising.onsager_internal_energy(β))
        @test Ising.onsager_heat_capacity(β) > 0
    end
end

@testset "Onsager internal energy and heat capacity at high and low temperature" begin
    # reference internal energy in high precision, using 2/π * K(m) = 1 / agm(1, sqrt(1 - m))
    function u_exact(β)
        setprecision(BigFloat, 256) do
            t = tanh(2big(β)); j = 2t^2 - 1
            a, g = big(1), abs(j)
            for _ in 1:30
                a, g = (a + g) / 2, sqrt(a * g)
            end
            Float64(-(1 + j / a) / t)
        end
    end
    for β in [10.0 .^ (-12:0.25:-1); 0.1:0.05:5]
        @test Ising.onsager_internal_energy(β) ≈ u_exact(β) rtol=1e-12
    end

    @test iszero(Ising.onsager_internal_energy(0))
    @test iszero(Ising.onsager_heat_capacity(0))
    @test Ising.onsager_internal_energy(Inf) == -2
    @test iszero(Ising.onsager_heat_capacity(Inf))
    # high temperature expansions, u = -2β + O(β^3), c = 2β^2 + O(β^4)
    for β in (1e-12, 1e-8, 1e-6)
        @test Ising.onsager_internal_energy(β) ≈ -2β rtol=1e-10
        @test Ising.onsager_heat_capacity(β) ≈ 2β^2 rtol=1e-10
    end
    # low temperature expansion, c = 64β^2 exp(-8β) (1 + O(exp(-4β)))
    for β in (10, 15, 20)
        @test Ising.onsager_heat_capacity(β) ≈ 64β^2 * exp(-8β) rtol=1e-12
    end
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
