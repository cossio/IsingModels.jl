bits2spins(spins::AbstractMatrix{Bool}) = Int8(2) * spins .- Int8(1)
spins2bits(spins::AbstractMatrix{Int8}) = BitMatrix(spins .> 0)

function neighbor_sum(spins::AbstractMatrix{Bool}, i::Int, j::Int)
    K, L = size(spins)
    @assert 1 ≤ i ≤ L && 1 ≤ j ≤ K
    top = spins[mod1(i + 1, K), j]
    bot = spins[mod1(i - 1, K), j]
    lef = spins[i, mod1(j + 1, L)]
    rig = spins[i, mod1(j - 1, L)]
    return 2 * (top + bot + lef + rig - 2)
end

function neighbor_sum_div_2(spins::AbstractMatrix{Bool}, i::Int, j::Int)
    K, L = size(spins)
    @assert 1 ≤ i ≤ L && 1 ≤ j ≤ K
    top = spins[mod1(i + 1, K), j]
    bot = spins[mod1(i - 1, K), j]
    lef = spins[i, mod1(j + 1, L)]
    rig = spins[i, mod1(j - 1, L)]
    return top + bot + lef + rig - 2
end

random_configuration(::Type{Bool}, L::Int, K::Int = L) = BitMatrix(rand(Bool, L, K))

function energy(spins::AbstractMatrix{Bool})
    E = 0
    for i = 1:N, j = 1:M
        E -= spins[i,j] * (spins[mod1(i + 1, N), j] + spins[i, mod1(j + 1, M)])
    end
    E = 4E + 4sum(spins) - length(spins)
    return E
end
