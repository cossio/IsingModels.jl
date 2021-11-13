#= Based on https://github.com/prcastro/IsingLite.jl/blob/master/src/wolff.jl =#

function recursion_wolff(
    spins::Matrix{Int8},   # Spin grid
    i::Int, j::Int,     # current position
    Padd::Real,
    cluster::BitMatrix = falses(size(spins)...) # cluster we are building
)
    cluster[i,j] = true
    for (x, y) in neighbors(i, j, size(spins)...)
        if cluster[x,y] == false && spins[x,y] == spins[i,j]
            if rand() < Padd
                recursion_wolff(spins, x, y, Padd, cluster)
            end
        end
    end
    return cluster
end

function wolff_cluster(spins::Matrix{Int8}, i::Int, j::Int, Padd::Real)
    cluster = falses(size(spins)...)
    cluster[i,j] = true
    queue = [(i,j)]
    while !isempty(queue)
        (i,j) = pop!(queue)
        for (x,y) in neighbors(i, j, size(spins)...)
            if cluster[x,y] == false && spins[x,y] == spins[i,j] && rand() < Padd
                cluster[x,y] = true
                push!(queue, (x,y))
            end
        end
    end
    return cluster
end

"""
    wolff!(spins, β, steps = 1; save_interval = 1)

Perfoms one or more Wolff MC steps from the configuration `spins`, at inverse
temperature `β`.
"""
function wolff!(
    spins::Matrix{Int8},
    β::Real, # inverse temperature
    steps::Integer = 1; # number of iterations
    save_interval::Int = 1
)
    m = zeros(steps)
    E = zeros(steps)
    m[1] = mean(spins)
    E[1] = energy(spins) / length(spins)
    spins_t = [copy(spins)]

    Padd = 1 - exp(-2β)
    for t in 2:steps
        i, j = rand.(Base.OneTo.(size(spins)))
        cluster = recursion_wolff(spins, i, j, Padd)

        # change in magnetization
        Δm = 2spins[i,j] * mean(cluster)
        m[t] = m[t - 1] - Δm

        # flip cluster
        spins .= (1 .- 2cluster) .* spins

        # compute new energy
        E[t] = energy(spins) / length(spins)

        if save_interval !== nothing && t % save_interval == 0
            push!(spins_t, copy(spins))
        end
    end
    return spins_t, m, E
end
