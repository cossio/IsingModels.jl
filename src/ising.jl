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
    energy(spins, h = 0)

Computes the energy of `spins` in the Ising 2-dimensional model.
"""
function energy(spins::AbstractMatrix{Int8})
    L, K = size(spins)
    E = 0
    for i = 1:L, j = 1:K
        E -= spins[i,j] * (spins[mod1(i + 1, L), j] + spins[i, mod1(j + 1, K)])
    end
    return E
end

energy(spins::AbstractMatrix{Int8}, h::Real) = energy(spins) - h * sum(spins)

"""
    neighbor_sum(spins, i, j)

Sum of spins in neighbor sites of (i, j).
"""
function neighbor_sum(spins::AbstractMatrix{Int8}, i::Int, j::Int)
    L, K = size(spins)
    @assert (1 ≤ i ≤ L) && (1 ≤ j ≤ K)
    top = spins[mod1(i + 1, L), j]
    bot = spins[mod1(i - 1, L), j]
    lef = spins[i, mod1(j + 1, K)]
    rig = spins[i, mod1(j - 1, K)]
    return top + bot + lef + rig
end

"""
    neighbor_sum_div_2(spins, i, j)

neighbor_sum(spins, i, j) is always even. Therefore here we compute the sum, divided by 2.
"""
function neighbor_sum_div_2(spins::AbstractMatrix{Int8}, i::Int, j::Int)
    L, K = size(spins)
    @assert (1 ≤ i ≤ L) && (1 ≤ j ≤ K)
    top = spins[mod1(i + 1, L), j] > 0
    bot = spins[mod1(i - 1, L), j] > 0
    lef = spins[i, mod1(j + 1, K)] > 0
    rig = spins[i, mod1(j - 1, K)] > 0
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
neighbors(spins::AbstractMatrix, i::Int, j::Int) = neighbors(size(spins)..., i, j)

"""
    random_configuration(L, K = L)

Generate a random spin configuration.
"""
random_configuration(::Type{Int8}, L::Int, K::Int = L) = rand((Int8(1), Int8(-1)), L, K)
random_configuration(L::Int, K::Int = L) = random_configuration(Int8, L, K)

"""
    random_magnetized_configuration(M, L, K = L)

Generate a random spin configuration with magnetization (sum over all spins) equal to `M`.
"""
function random_magnetized_configuration(M::Int, L::Int, K::Int = L)
    N = L * K
    @assert -N ≤ M ≤ N
    @assert iseven(N + M) # 2 * number of +1 spins
    @assert iseven(N - M) # 2 * number of -1 spins
    Np = (N + M) ÷ 2 # number of +1 spins
    Nm = (N - M) ÷ 2 # number of +1 spins
    @assert Np + Nm == N
    spins = zeros(Int8, L, K)
    p = randperm(L * K)
    spins[p[1:Np]] .= +1
    spins[p[(end - Nm + 1):end]] .= -1
    @assert all((spins .== 1) .| (spins .== -1))
    return spins
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
distance_matrix(spins::AbstractMatrix) = distance_matrix(size(spins)...)
