function wolff_cluster(s::IsingArray{D}, i::CartesianIndex{D}, Padd::Real = 1) where {D}
    cluster = falses(size(s)...)
    cluster[i] = true
    queue = [i]
    while !isempty(queue)
        j = pop!(queue)
        for x in neighbors(CartesianIndices(s), j)
            if !cluster[x] && s[x] == s[j] && rand() < Padd
                cluster[x] = true
                push!(queue, x)
            end
        end
    end
    return cluster
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
    for t in 2:steps
        wolff_step!(σ; t = t, M = M, E = E, Padd = Padd)
        if t ∈ 1:save_interval:steps
            selectdim(σ_t, ndims(σ) + 1, cld(t, save_interval)) .= σ
        end
    end
    return σ_t, M, E
end

function wolff_step!(σ::BitArray; t, M, E, Padd)
    i = rand(CartesianIndices(σ))
    cluster = wolff_cluster(σ, i, Padd)
    cluster_size = sum(cluster)
    ΔM = -2spin(σ[i]) * cluster_size
    flip_cluster!(σ, cluster)
    M[t] = M[t - 1] + ΔM
    E[t] = energy(σ)
    return cluster_size
end

function flip_cluster!(σ::BitArray, cluster::BitArray)
    σ .= ifelse.(cluster, flip.(σ), σ)
    return σ
end

wolff_padd(β::Real) = -expm1(-2β)
