function wolff_cluster(s::AbstractMatrix, i::Int, j::Int, Padd::Real = 1)
    cluster = falses(size(s)...)
    cluster[i,j] = true
    queue = [(i,j)]
    while !isempty(queue)
        (i,j) = pop!(queue)
        for (x,y) in neighbors(s, i, j)
            if !cluster[x,y] && s[x,y] == s[i,j] && rand() < Padd
                cluster[x,y] = true
                push!(queue, (x,y))
            end
        end
    end
    return cluster
end

"""
    wolff!(spins, β; steps = 1, save_interval = 1)

Perfoms one or more Wolff MC steps from the configuration `spins`, at inverse
temperature `β`.
"""
function wolff!(
    spins::AbstractMatrix,
    β::Real; # inverse temperature
    steps::Integer = 1, # number of iterations
    save_interval::Int = 1
)
    M = zeros(Int, steps)
    E = zeros(Int, steps)
    M[1] = magnetization(spins)
    E[1] = energy(spins)
    spins_t = zeros(eltype(spins), size(spins)..., length(1:save_interval:steps))
    spins_t[:,:,1] .= spins

    Padd = wolff_padd(β)
    for t in 2:steps
        wolff_step!(spins; t = t, M = M, E = E, Padd = Padd)
        if t ∈ 1:save_interval:steps
            spins_t[:, :, cld(t, save_interval)] .= spins
        end
    end
    return spins_t, M, E
end

function wolff_step!(s::AbstractMatrix; t, M, E, Padd)
    i, j = rand.(Base.OneTo.(size(s)))
    cluster = wolff_cluster(s, i, j, Padd)
    cluster_size = sum(cluster)
    M[t] = M[t - 1] - 2ising(s[i,j]) * cluster_size
    flip_cluster!(s, cluster)
    E[t] = energy(s)
    return cluster_size
end

function flip_cluster!(s::AbstractMatrix{Int8}, cluster::AbstractMatrix{Bool})
    s .= (1 .- 2cluster) .* s
    return s
end

function flip_cluster!(σ::AbstractMatrix{Bool}, cluster::AbstractMatrix{Bool})
    σ .= ifelse.(cluster, flip.(σ), σ)
    return σ
end

wolff_padd(β::Real) = -expm1(-2β)
