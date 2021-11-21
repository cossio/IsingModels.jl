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
    energy(spins)

Computes the energy of `spins` in the Ising 2-dimensional model.
"""
function energy(spins::AbstractMatrix{Int8})
    N, M = size(spins)
    E = 0
    for i = 1:N, j = 1:M
        E -= spins[i,j] * (spins[mod1(i + 1, N), j] + spins[i, mod1(j + 1, M)])
    end
    return E
end

"""
    neighbor_sum(spins, i, j)

Sum of spins in neighbor sites of (i, j).
"""
function neighbor_sum(spins::AbstractMatrix{Int8}, i::Int, j::Int)
    S = 0
    for d in (-1, 1)
        S += spins[mod1(i + d, size(spins, 1)), j]
        S += spins[i, mod1(j + d, size(spins, 2))]
    end
    return S
end

"""
    neighbors(i, j, L, K = L)

Returns the list of neighbors of site `(i,j)` in the 2-dimensional lattice grid with sides
`L x K`. That is, the sites `((i+1,j), (i-1,j), (i,j+1), (i,j-1))`, but considering the
periodic bounary conditions at the edges.
"""
function neighbors(i::Int, j::Int, L::Int, K::Int = L)
    @assert 1 ≤ i ≤ L && 1 ≤ j ≤ K
    return ((mod1(i + 1, L), j),
            (mod1(i - 1, L), j),
            (i, mod1(j + 1, K)),
            (i, mod1(j - 1, K)))
end

"""
    random_configuration(L, K = L)

Generate a random spin configuration.
"""
random_configuration(L::Int, K::Int = L) = rand((Int8(1), Int8(-1)), L, K)

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

"""
    kramers_wannier(β)

Returns the Kramers-Wannier dual inverse-temperature of β.
"""
kramers_wannier(β::Real) = -log(tanh(β)) / 2
