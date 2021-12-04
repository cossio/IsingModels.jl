"""
    adjacency_matrix(L, K = L)

Adjacency matrix of a L x K rectangular grid lattice.
"""
function adjacency_matrix(L::Int, K::Int = L)
    A = falses(L, K, L, K)
    for i in 1:L, j in 1:K, d in (-1, 1)
        A[i, j, mod1(i + d, L), j] = true
        A[i, j, i, mod1(j + d, K)] = true
    end
    return reshape(A, L * K, L * K)
end

"""
    ising(s)

Returns `+1` or `-1`, from different representations of a spin.
"""
ising(s::Signed) = ifelse(s > 0, oftype(s, 1), oftype(s, -1))
ising(σ::Bool) = Int8(2σ - 1)

"""
    binary(s)

Returns `true` or `false`, from different representations of a spin.
"""
binary(s::Signed) = s > 0
binary(σ::Bool) = σ

"""
    energy(spins, h = 0)

Computes the energy of `spins` in the Ising 2-dimensional model.
"""
function energy(s::AbstractMatrix{T}) where {T<:Signed}
    L, K = size(s)
    τ = ising.(s)
    return -sum(τ .* (τ[mod1.(2:(L + 1), L), :] .+ τ[:, mod1.(2:(K + 1), K)]))
end

energy(s::AbstractMatrix{T}, h::Real) where {T<:Signed} = energy(s) - h * sum(ising, s)

function energy(σ::AbstractMatrix{Bool})
    L, K = size(σ)
    couplings = σ .* (σ[mod1.(2:(L + 1), L), :] .+ σ[:, mod1.(2:(K + 1), K)])
    return -4sum(couplings) + 8sum(σ) - 2length(σ)
end

energy(σ::AbstractMatrix{Bool}, h::Real) = energy(σ) - 2h * sum(σ) + h * length(σ)

magnetization(s::AbstractMatrix{T})  where {T<:Signed} = sum(ising, s)
magnetization(σ::AbstractMatrix{Bool}) = 2sum(σ) - length(σ)

flip(s::Signed) = ifelse(s > 0, oftype(s, -1), oftype(s, 1))
flip(σ::Bool) = !σ

"""
    neighbor_sum(spins, i, j)

Sum of spins in neighbor sites of (i, j).
"""
function neighbor_sum(s::AbstractMatrix, i::Int, j::Int)
    L, K = size(s)
    @assert (1 ≤ i ≤ L) && (1 ≤ j ≤ K)
    top = ising(s[mod1(i + 1, L), j])
    bot = ising(s[mod1(i - 1, L), j])
    lef = ising(s[i, mod1(j + 1, K)])
    rig = ising(s[i, mod1(j - 1, K)])
    return top + bot + lef + rig
end

"""
    neighbor_sum_div_2(spins, i, j)

neighbor_sum(spins, i, j) is always even. Therefore here we compute the sum, divided by 2.
"""
function neighbor_sum_div_2(s::AbstractMatrix, i::Int, j::Int)
    L, K = size(s)
    @assert (1 ≤ i ≤ L) && (1 ≤ j ≤ K)
    top = s[mod1(i + 1, L), j] > 0
    bot = s[mod1(i - 1, L), j] > 0
    lef = s[i, mod1(j + 1, K)] > 0
    rig = s[i, mod1(j - 1, K)] > 0
    return top + bot + lef + rig - 2
end

"""
    neighbors(L, K [= L], i, j)
    neighbors(spins, i, j)

Returns the list of neighbors of site `(i,j)` in the 2-dimensional lattice grid with sides
`L x K`. That is, the sites `((i+1,j), (i-1,j), (i,j+1), (i,j-1))`, but considering the
periodic bounary conditions at the edges.
"""
function neighbors(L::Int, K::Int, i::Int, j::Int)
    @assert 1 ≤ i ≤ L && 1 ≤ j ≤ K
    return ((mod1(i + 1, L), j),
            (mod1(i - 1, L), j),
            (i, mod1(j + 1, K)),
            (i, mod1(j - 1, K)))
end

neighbors(L::Int, i::Int, j::Int) = neighbors(L, L, i, j)
neighbors(s::AbstractMatrix, i::Int, j::Int) = neighbors(size(s)..., i, j)

"""
    random_configuration(L, K = L)

Generate a random spin configuration.
"""
function random_configuration(::Type{T}, L::Int, K::Int = L) where {T<:Signed}
    return rand(T.((-1, 1)), L, K)
end
random_configuration(::Type{Bool}, L::Int, K::Int = L) = bitrand(L, K)
random_configuration(L::Int, K::Int = L) = random_configuration(Int8, L, K)

"""
    random_magnetized_configuration(M, L, K = L)

Generate a random spin configuration with magnetization (sum over all spins) equal to `M`.
You must make sure that the given value `M` is feasible.
"""
function random_magnetized_configuration(::Type{Bool}, M::Int, L::Int, K::Int = L)
    N = L * K
    @assert -N ≤ M ≤ N
    @assert iseven(N + M) # 2 * number of +1 spins
    @assert iseven(N - M) # 2 * number of -1 spins
    Np = (N + M) ÷ 2 # number of +1 spins
    Nm = (N - M) ÷ 2 # number of +1 spins
    @assert Np + Nm == N
    spins = falses(L, K)
    p = randperm(L * K)
    spins[p[1:Np]] .= true
    spins[p[(end - Nm + 1):end]] .= false
    @assert all((spins .== true) .| (spins .== false))
    return spins
end

function random_magnetized_configuration(
    ::Type{T}, M::Int, L::Int, K::Int = L
) where {T<:Signed}
    σ = random_magnetized_configuration(Bool, M, L, K)
    return T.(ising.(σ))
end

function random_magnetized_configuration(M::Int, L::Int, K::Int = L)
    return random_magnetized_configuration(Int8, M, L, K)
end

"""
    distance_tensor(L, K = L)

Returns a `LxKxLxK` tensor `dist`, such that the entry
`dist[i1,j1,i2,j2]` gives the distance between sites
`(i1,j1)` and `(i2,j2)` in the lattice.
"""
function distance_tensor(L::Int, K::Int = L)
    d1 = (reshape(1:L, L, 1) .- reshape(1:L, 1, L))
    d2 = (reshape(1:K, K, 1) .- reshape(1:K, 1, K))
    @. d1 = min(d1^2, (d1 - L)^2, (d1 + L)^2)
    @. d2 = min(d2^2, (d2 - K)^2, (d2 + K)^2)
    D = reshape(d1,L,1,L,1) .+ reshape(d2,1,K,1,K)
    return sqrt.(D)
end

"""
    distance_matrix(L, K = L)
    distance_matrix(spins)

Returns a `LxK` matrix `dist`, such that the entry
`dist[i,j]` gives the distance between sites
`(0,0)` and `(i,j)` in the lattice.
Takes into account the periodicity of the lattice.
"""
function distance_matrix(L::Int, K::Int = L)
    x = 0:(L - 1)
    y = 0:(K - 1)
    dx = @. min(x^2, (x - L)^2, (x + L)^2)
    dy = @. min(y^2, (y - K)^2, (y + K)^2)
    d = reshape(dx, L, 1) .+ reshape(dy, 1, K)
    return sqrt.(d)
end

distance_matrix(s::AbstractMatrix) = distance_matrix(size(s)...)
