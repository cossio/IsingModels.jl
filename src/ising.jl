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
function neighbor_sum(spins::Matrix{Int8}, i::Int, j::Int)
    S = 0
    for d in (-1, 1)
        S += spins[mod1(i + d, size(spins, 1)), j]
        S += spins[i, mod1(j + d, size(spins, 2))]
    end
    return S
end

"""
    neighbors(i, j, N, M = N)

Returns the list of neighbors of site `(i,j)` in the 2-dimensional lattice grid with sides
`N x M`. That is, the sites `((i+1,j), (i-1,j), (i,j+1), (i,j-1))`, but considering the
periodic bounary conditions at the edges.
"""
function neighbors(i::Int, j::Int, N::Int, M::Int = N)
    @assert 1 ≤ i ≤ N && 1 ≤ j ≤ M
    return ((mod1(i + 1, N), j),
            (mod1(i - 1, N), j),
            (i, mod1(j + 1, M)),
            (i, mod1(j - 1, M)))
end

"""
    random_configuration(L, K = L)

Generate a random spin configuration.
"""
random_configuration(L::Int, K::Int = L) = rand((Int8(1), Int8(-1)), L, K)
