function wolff_cluster(s::IsingArray{D}, i::CartesianIndex{D}, Padd::Real = 1) where {D}
    cluster = falses(size(s)...)
    wolff_cluster!(cluster, s, i, Padd)
    return cluster
end

# Grows the cluster from the seed `i`, marking its sites in `cluster`, which must be
# all `false` on entry. Returns the list of cluster sites.
function wolff_cluster!(
    cluster::IsingArray{D}, s::IsingArray{D}, i::CartesianIndex{D}, Padd::Real = 1
) where {D}
    cluster[i] = true
    sites = [i]
    queue = [i]
    while !isempty(queue)
        j = pop!(queue)
        for x in neighbors(CartesianIndices(s), j)
            if !cluster[x] && s[x] == s[j] && rand() < Padd
                cluster[x] = true
                push!(queue, x)
                push!(sites, x)
            end
        end
    end
    return sites
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
    sites = wolff_cluster!(cluster, σ, i, Padd)
    cluster_size = length(sites)
    ΔM = -2spin(σ[i]) * cluster_size
    ΔE = flip_cluster_ΔE(σ, cluster, sites)
    flip_cluster!(σ, sites)
    cluster[sites] .= false
    M[t] = M[t - 1] + ΔM
    E[t] = E[t - 1] + ΔE
    return cluster_size
end

function flip_cluster!(σ::BitArray, cluster::BitArray)
    σ .= ifelse.(cluster, flip.(σ), σ)
    return σ
end

function flip_cluster!(σ::BitArray, sites::AbstractVector{<:CartesianIndex})
    for i in sites
        σ[i] = flip(σ[i])
    end
    return σ
end

# Energy change from flipping the cluster, computed before the flip. Only the bonds
# between a cluster site and an outside neighbor change sign. Like `energy`, this
# counts a bond twice along a dimension of size 2, where both neighbors coincide.
function flip_cluster_ΔE(σ::IsingArray{D}, cluster::IsingArray{D}, sites) where {D}
    ΔE = 0
    for i in sites, j in neighbors(CartesianIndices(σ), i)
        if !cluster[j]
            ΔE += 2spin(σ[i]) * spin(σ[j])
        end
    end
    return ΔE
end

wolff_padd(β::Real) = -expm1(-2β)
