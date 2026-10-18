function wolff_cluster(s::IsingArray{D}, i::CartesianIndex{D}, Padd::Real = 1) where {D}
    σ = copy(s)
    wolff_flip!(falses(size(s)), σ, i, Padd)
    return σ .!= s # the flipped sites
end

# Grows a Wolff cluster from the seed `i` and flips it. Returns the cluster size and the
# energy change. `queued` is a work buffer of the size of `σ`, which must be all `false`:
# it marks the cluster sites whose neighbors are not scanned yet, and is left all `false`.
function wolff_flip!(
    queued::IsingArray{D}, σ::IsingArray{D}, i::CartesianIndex{D}, Padd::Real = 1
) where {D}
    @assert size(queued) == size(σ)
    σi = σ[i] # spin of the cluster sites before the flip
    queued[i] = true
    queue = [i]
    cluster_size = 0
    up = 0 # up neighbors of the cluster sites, each counted just before the site flips
    while !isempty(queue)
        j = pop!(queue)
        for x in neighbors(CartesianIndices(σ), j)
            σx = @inbounds σ[x]
            up += σx
            # sites already in the cluster are either queued, or flipped to !σi
            if !(@inbounds queued[x]) && σx == σi && rand() < Padd
                @inbounds queued[x] = true
                push!(queue, x)
            end
        end
        @inbounds σ[j] = flip(σ[j])
        @inbounds queued[j] = false
        cluster_size += 1
    end
    # This flips the cluster one site at a time. Flipping j changes the energy by
    # 2 s_j Σ_x s_x over its 2D neighbors x, where s_j = spin(σi), and Σ_x s_x is the
    # number of up neighbors minus down neighbors at that time. Along a dimension of size 1,
    # j is its own neighbor, but this bond doesn't change, so its 4 is subtracted.
    n1 = count(==(1), size(σ))
    ΔE = 2spin(σi) * (2up - 2D * cluster_size) - 4n1 * cluster_size
    return cluster_size, ΔE
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
    queued = falses(size(σ)) # work buffer, reused across steps
    for t in 2:steps
        wolff_step!(σ; t = t, M = M, E = E, Padd = Padd, queued = queued)
        if t ∈ 1:save_interval:steps
            selectdim(σ_t, ndims(σ) + 1, cld(t, save_interval)) .= σ
        end
    end
    return σ_t, M, E
end

# `queued` is the work buffer of `wolff_flip!`, which can be reused across steps.
function wolff_step!(σ::IsingArray; t, M, E, Padd, queued::IsingArray = falses(size(σ)))
    i = rand(CartesianIndices(σ))
    ΔM = -2spin(σ[i])
    cluster_size, ΔE = wolff_flip!(queued, σ, i, Padd)
    M[t] = M[t - 1] + ΔM * cluster_size
    E[t] = E[t - 1] + ΔE
    return cluster_size
end

function flip_cluster!(σ::IsingArray, cluster::IsingArray)
    σ .= ifelse.(cluster, flip.(σ), σ)
    return σ
end

wolff_padd(β::Real) = -expm1(-2β)
