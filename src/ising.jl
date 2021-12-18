const IsingIndices{D} = CartesianIndices{D, NTuple{D, Base.OneTo{Int}}}
const IsingArray{D} = AbstractArray{Bool,D}
const IsingVector = IsingArray{1}

"""
    spin(σ)

Returns `+1` or `-1` spin value, from a boolean `σ` representation.
"""
spin(σ::Bool) = 2σ - 1

"""
    binary(s)

Converts to Boolean representation.
"""
binary(s::Signed) = s > 0

"""
    flip(σ)

Flips a spin.
"""
flip(σ::Bool) = !σ

"""
    magnetization(σ)

Computes the magnetization of a set of spins.
"""
magnetization(σ::IsingArray) = 2sum(σ) - length(σ)

"""
    energy(σ, h = 0; f)

Computes the energy of `σ` in the Ising model.
"""
function energy(σ::IsingArray, h::Real = 0; f = nothing)
    sz = CartesianIndices(σ)
    N = length(σ)
    c = sum(σ[i] * σ[j] for i in sz for j in neighbors_p(sz, i))
    M = magnetization(σ)
    E = -4c + (2M + N) * ndims(σ) - h * M
    if isnothing(f)
        return E
    else
        return E + f(M)
    end
end

"""
    sum_neighbors(σ, i)

Sum of spins in neighbor sites of i.
"""
function sum_neighbors(σ::IsingArray{D}, i::CartesianIndex{D}) where {D}
    @boundscheck @assert i ∈ CartesianIndices(σ)
    c = sum(map(j -> σ[j], neighbors(CartesianIndices(σ), i)))
    return 2 * (c - ndims(σ))
end

function neighbor(sz::IsingIndices{D}, i::CartesianIndex{D}, dim::Int, Δ::Int) where {D}
    return CartesianIndex(ntuple(d -> d == dim ? mod1(i[d] + Δ, size(sz, d)) : i[d], D))
end

function neighbors(sz::IsingIndices{D}, i::CartesianIndex{D}) where {D}
    return (neighbors_p(sz, i)..., neighbors_m(sz, i)...)
end

function neighbors_p(sz::CartesianIndices{D}, i::CartesianIndex{D}) where {D}
    return ntuple(d -> neighbor(sz, i, d, 1), ndims(sz))
end

function neighbors_m(sz::CartesianIndices, i::CartesianIndex)
    return ntuple(d -> neighbor(sz, i, d, -1), ndims(sz))
end

"""
    random_magnetized_configuration!(σ, M)

Generate a random spin configuration with magnetization equal to `M`.
You must make sure that the given value `M` is feasible.
"""
function random_magnetized_configuration!(σ::IsingArray, M::Int)
    N = length(σ)
    @assert -N ≤ M ≤ N
    @assert iseven(N + M) # 2 * number of +1 spins
    P = (N + M) ÷ 2 # number of +1 spins
    r = randperm(N)
    σ[r[1:P]] .= true
    σ[r[(P + 1):end]] .= false
    return σ
end

function ising_distance(
    sz::IsingIndices{D}, i::CartesianIndex{D}, j::CartesianIndex{D}
) where {D}
    @boundscheck @assert i ∈ sz && j ∈ sz
    ds = ising_axis_distance.(size(sz), Tuple(i - j))
    return norm(ds)
end

function ising_axis_distance(L::Int, i::Int, j::Int)
    d = abs(i - j)
    if 2d ≤ L
        return d
    else
        return L - d
    end
end
