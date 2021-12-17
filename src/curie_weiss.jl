curie_weiss_magnetization(spins::AbstractVector) = sum(ising, spins)

function curie_weiss_energy(spins::AbstractVector, h::Real = false; f = nothing)
    N = length(spins)
    M = curie_weiss_magnetization(spins)
    E = 1/2 * (1 - M^2 / N) - h * M
    if isnothing(f)
        return E
    else
        return E + f(M)
    end
end

function curie_weiss_random_configuration!(spins::AbstractVector{Int8}, M::Int)
    N = length(spins)
    @assert iseven(N + M)
    Np = (N + M) ÷ 2
    p = randperm(N)
    spins[p[begin:Np]] .= 1
    spins[p[(Np + 1):end]] .= -1
    @assert sum(spins) == M
    return spins
end

function curie_weiss_random_configuration(
    N::Int, β::Real, h::Real = false; f = nothing, B::Int = 1
)
    Ms = curie_weiss_magnetization_rand(N, β, h; f = f, B = B)
    @assert all(iseven, M .+ N)
    spins = zeros(Int8, N, B)
    for (b, M) in zip(1:B, Ms)
        curie_weiss_random_configuration!(view(spins, :, b), M)
    end
    return spins
end

@doc raw"""
    curie_weiss_random_magnetizations(N, β, h = 0; f = nothing, B = 1)

Generates a random magnetization of the Curie-Weiss model, that is, a value of

```math
M = \sum_{i=1}^N s_i
```

that follows the statistics of the Curie-Weiss model, with energy function:

```math
E = -\sum_{i < j} s_i s_j - h * \sum_i s_i + f(M)
```

By default f(M) = 0.

`B` controls the number of samples generated.
"""
function curie_weiss_random_magnetizations(
    N::Int, β::Real, h::Real = false; f = nothing, B::Int = 1
)
    E(M) = SpecialFunctions.logabsbinomial(N, M) + M^2 / 2N + h * M - f(M)
    p = LogExpFunctions.softmax(-β * E.(-N:2:N))
    return 2rand(Distributions.Categorical(p), B) .- N
end
