"""
    mf_energy(σ, h = 0; f = nothing)
"""
function mf_energy(σ::IsingVector, h::Real = false; f = nothing)
    N = length(σ)
    M = magnetization(σ)
    E = 1/2 * (1 - M^2 / N) - h * M
    if isnothing(f)
        return E
    else
        return E + f(M)
    end
end

"""
    mf_random_configurations(N, β, h = 0; f, B = 1)

Generates `B` random configurations of the Curie-Weiss model.
"""
function mf_random_configurations(
    N::Int, β::Real, h::Real = false; f = nothing, B::Int = 1
)
    Ms = mf_random_magnetizations(N, β, h; f = f, B = B)
    @assert all(iseven, M .+ N)
    σ = falses(N, B)
    for (b, M) in zip(1:B, Ms)
        mf_random_configuration!(view(σ, :, b), M)
    end
    return σ
end

"""
    mf_random_magnetized_configuration!(σ, M)

Random configuration with magnetization equal to `M`.
"""
function mf_random_magnetized_configuration!(σ::IsingArray{1}, M::Int)
    N = length(σ)
    @assert iseven(N + M)
    Np = (N + M) ÷ 2

    p = randperm(N)
    σ[p[1:Np]] .= true
    σ[p[(Np + 1):end]] .= false

    return σ
end

"""
    mf_random_magnetized_configuration(N, M)

Random Curie-Weiss configuration with `N` spind, with magnetization equal to `M`.
"""
function mf_random_magnetized_configuration(N::Int, M::Int)
    @assert iseven(N + M)
    σ = BitVector(undef, N)
    mf_random_magnetized_configuration!(σ, M)
    return σ
end

@doc raw"""
    mf_random_magnetizations(N, β, h = 0; f = nothing, B = 1)

Generates a random magnetization of the Curie-Weiss model, that is, a value of

```math
M = \sum_{i=1}^N s_i
```

that follows the statistics of the mean-field (MF) Curie-Weiss model, with energy function:

```math
E = -\sum_{i < j} s_i s_j - h * \sum_i s_i + f(M)
```

By default f(M) = 0.

`B` controls the number of samples generated.
"""
function mf_random_magnetizations(
    N::Int, β::Real, h::Real = false; f = nothing, B::Int = 1
)
    E(M) = SpecialFunctions.logabsbinomial(N, M) + M^2 / 2N + h * M - f(M)
    p = LogExpFunctions.softmax(-β * E.(-N:2:N))
    Np = rand(Distributions.Categorical(p), B)
    return 2Np .- N
end
