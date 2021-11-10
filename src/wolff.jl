#= Based on https://github.com/prcastro/IsingLite.jl/blob/master/src/wolff.jl =#

function recursion_wolff(
    spins::Matrix{Int8},   # Spin grid
    Padd::Real,
    i::Int, j::Int,     # current position
    cluster::BitMatrix = falses(size(spin)...) # cluster we are building
)
    cluster[i,j] = true
    for (x, y) in neighbors(i, j, size(spins)...)
        if cluster[x,y] == false && spins[x,y] == spins[i,j]
            if rand() < Padd
                recursion_wolff(spins, Padd, x, y, cluster)
            end
        end
    end
    return cluster
end

function wolff!(
    spins::Matrix{Int8},
    β::Real;  # inverse temperature
    iters::Integer = 30,   # Number of iterations
)
    N, M = size(spins)
    Padd = 1 - exp(-2β)
    m = Float64[]
    for t in 1:iters
        i, j = rand(1:N), rand(1:M)
        cluster = recursion_wolff(spins, i, j, Padd)
        ΔE = -2h * dot(spins, cluster)
        # Change spin accordingly
	    if ΔE ≤ 0 || rand() < exp(-ΔE)
            flip!(spins, cluster)
        end
        Δm =
        push!(m, mean(spins))
    end

    return m
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

function flip!(grid::Array{Int, 2}, cluster::BitArray{2})
    for j in 1:size(grid, 2), i in 1:size(grid, 1)
        if cluster[i, j] == true
            grid[i, j] *= -1
        end
    end
end
