include("init.jl")

@test Ising.ising(true) == 1
@test Ising.ising(false) == -1
@test Ising.ising(Int8(+1)) == Int8(+1)
@test Ising.ising(Int8(-1)) == Int8(-1)

for σ in (true, false)
    @test Ising.binary(Ising.ising(σ)) == σ
end

for s in Int8[-1, 1]
    @test Ising.ising(Ising.binary(s)) == s
end

@test Ising.binary(true) == true
@test Ising.binary(false) == false
@test Ising.binary(Int8(+1)) == true
@test Ising.binary(Int8(-1)) == false

σ = bitrand(5, 7)
@test Ising.binary.(Ising.ising.(σ)) == σ
@test Ising.energy(σ) == Ising.energy(Ising.ising.(σ))
@test Ising.magnetization(σ) == Ising.magnetization(Ising.ising.(σ))

s = rand(Int8.((-1, 1)), 5, 7)
@test Ising.ising.(Ising.binary.(s)) == s
@test Ising.energy(s) == Ising.energy(Ising.binary.(s))
@test Ising.magnetization(s) == Ising.magnetization(Ising.binary.(s))

h = randn()
@test Ising.energy(s, h) ≈ Ising.energy(Ising.binary.(s), h)
@test Ising.energy(σ, h) ≈ Ising.energy(Ising.ising.(σ), h)
