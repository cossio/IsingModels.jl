function wolff_cluster(s::IsingArray{D}, i::CartesianIndex{D}, Padd::Real = 1) where {D}
    cluster = falses(size(s)...)
    wolff_flip!(cluster, copy(s), i, Padd) # flips a copy, leaving `s` unchanged
    return cluster
end

# Grows a Wolff cluster from the seed `i` and flips it, marking its sites in `cluster`,
# which must be all `false` on entry. Returns the list of cluster sites and the energy
# change of the flip.
function wolff_flip!(
    cluster::IsingArray{D}, σ::IsingArray{D}, i::CartesianIndex{D}, Padd::Real = 1
) where {D}
    σi = σ[i] # spin of all cluster sites before the flip
    # The flip changes the sign of the bonds between the cluster and outside sites, so
    # ΔE = 2(n_same - n_opp), counting these bonds by whether the outside spin is the same
    # or opposite. Over the 2D bonds of each cluster site, n_same = 2D|C| - n_opp - n_in,
    # where n_in counts bonds with both ends in the cluster, once from each end. Like
    # `energy`, this counts a bond twice along a dimension of size 2.
    n_opp = n_in = 0
    cluster[i] = true
    sites = [i]
    queue = [i]
    while !isempty(queue)
        j = pop!(queue)
        for x in neighbors(CartesianIndices(σ), j)
            if cluster[x]
                if σ[x] != σi # x was processed and flipped before j, so count both ends now
                    n_in += 2
                elseif x == j # j is its own neighbor along a dimension of size 1
                    n_in += 1
                end
            elseif σ[x] != σi
                n_opp += 1
            elseif rand() < Padd
                cluster[x] = true
                push!(queue, x)
                push!(sites, x)
            end
        end
        σ[j] = flip(σ[j])
    end
    n_same = 2D * length(sites) - n_opp - n_in
    return sites, 2(n_same - n_opp)
end

"""
    wolff!(σ, β; steps = 1, save_interval = 1)

Perfoms one or more Wolff MC steps from the configuration `σ`, at inverse
temperature `β`.
"""
function wolff!(
    σ::IsingArray,
    β::Real; # inverse temperature
    steps::Integer = 1, # number of iterations
    save_interval::Int = 1
)
    M = zeros(Int, steps)
    E = zeros(Int, steps)
    M[1] = magnetization(σ)
    E[1] = energy(σ)
    σ_t = falses(size(σ)..., length(1:save_interval:steps))
    selectdim(σ_t, ndims(σ) + 1, 1) .= σ

    Padd = wolff_padd(β)
    cluster = falses(size(σ)) # reused across steps
    for t in 2:steps
        wolff_step!(σ; t = t, M = M, E = E, Padd = Padd, cluster = cluster)
        if t ∈ 1:save_interval:steps
            selectdim(σ_t, ndims(σ) + 1, cld(t, save_interval)) .= σ
        end
    end
    return σ_t, M, E
end

# `cluster` is a work buffer of the size of `σ`. It must be all `false` on entry, and is
# left all `false` on exit, so it can be reused across steps.
function wolff_step!(σ::BitArray; t, M, E, Padd, cluster::BitArray = falses(size(σ)))
    i = rand(CartesianIndices(σ))
    ΔM = -2spin(σ[i])
    sites, ΔE = wolff_flip!(cluster, σ, i, Padd)
    cluster[sites] .= false
    cluster_size = length(sites)
    M[t] = M[t - 1] + ΔM * cluster_size
    E[t] = E[t - 1] + ΔE
    return cluster_size
end

function flip_cluster!(σ::BitArray, cluster::BitArray)
    σ .= ifelse.(cluster, flip.(σ), σ)
    return σ
end

wolff_padd(β::Real) = -expm1(-2β)
