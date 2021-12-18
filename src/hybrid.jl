"""
    hybrid!(σ, β; steps = 1, save_interval = length(σ), local_steps = length(σ))

Hybrid sampler, performing `local_steps` of Metropolis sampling, then one Wolff cluster
move, then another `local_steps` of Metropolis sampling, one more Wolff cluster move,
and so on.
"""
function hybrid!(
    σ::IsingArray,
    β::Real;
    steps::Int = 1,
    save_interval::Int = length(σ),
    local_steps::Int = length(σ) # Metropolis step per Wolff step
)
    @assert steps ≥ 1
    @assert save_interval ≥ 1
    @assert local_steps ≥ 1

    #= Track history of magnetization and energy =#
    M = zeros(Int, steps)
    E = zeros(Int, steps)

    M[1] = magnetization(σ)
    E[1] = energy(σ)

    #= Track the history of configurations only every 'save_interval' steps. =#
    σ_t = falses(size(σ)..., length(1:save_interval:steps))
    σ_t[:,:,1] .= σ

    Padd = wolff_padd(β)

    for t ∈ 2:steps
        if t ∈ 1:local_steps:steps # Wolff step
            wolff_step!(σ; t = t, M = M, E = E, Padd = Padd)
        else # Metropolis step
            metropolis_step!(σ; t = t, M = M, E = E, β = β)
        end
        if t ∈ 1:save_interval:steps
            σ_t[:, :, cld(t, save_interval)] .= σ
        end
    end

    return σ_t, M, E
end

mutable struct HybridStats
    local_flip::Float64
    wolff_flip::Float64
    local_time::Float64
    wolff_time::Float64
end
HybridStats() = HybridStats(0, 0, 0, 0)

"""
    dynamic_hybrid!(σ, β; steps, save_interval)

Same as `hybrid!`, but adjusts numbers of Metropolis and Wolff steps dynamically.
"""
function dynamic_hybrid!(
    σ::IsingArray,
    β::Real;
    steps::Int = 1,
    save_interval::Int = length(σ),
    hybrid_stats::HybridStats = HybridStats()
)
    @assert steps ≥ 1
    @assert save_interval ≥ 1

    #= Track history of magnetization and energy =#
    M = zeros(Int, steps)
    E = zeros(Int, steps)

    M[1] = magnetization(σ) # magnetization
    E[1] = energy(σ)

    #= Track the history of configurations only every 'save_interval' steps. =#
    σ_t = falses(size(σ)..., length(1:save_interval:steps))
    σ_t[:,:,1] .= σ

    Padd = wolff_padd(β)

    for t ∈ 2:steps
        if hybrid_decide(hybrid_stats, length(σ))
            hybrid_stats.wolff_time += @elapsed begin
                flipped = wolff_step!(σ; t = t, M = M, E = E, Padd = Padd)
                hybrid_stats.wolff_flip += flipped
            end
        else
            hybrid_stats.local_time += @elapsed begin
                flipped = metropolis_step!(σ; t = t, M = M, E = E, β = β)
                hybrid_stats.local_flip += flipped
            end
        end
        if t ∈ 1:save_interval:steps
            σ_t[:, :, cld(t, save_interval)] .= σ
        end
    end

    return σ_t, M, E
end

function hybrid_decide(hybrid_stats::HybridStats, N::Integer)
    # return true -> do Wolff, else do Metropolis
    DO_WOLFF = true
    DO_LOCAL = false

    if iszero(hybrid_stats.wolff_time) || iszero(hybrid_stats.wolff_flip)
        return DO_WOLFF
    elseif iszero(hybrid_stats.local_time) || iszero(hybrid_stats.local_flip)
        return DO_LOCAL
    elseif hybrid_stats.wolff_flip * N < hybrid_stats.local_flip
        return DO_WOLFF
    elseif hybrid_stats.local_flip * N < hybrid_stats.wolff_flip
        return DO_LOCAL
    end

    wolff_rate = hybrid_stats.wolff_flip / hybrid_stats.wolff_time
    local_rate = hybrid_stats.local_flip / hybrid_stats.local_time

    Pwolff = wolff_rate / (wolff_rate + local_rate)
    if rand() < Pwolff
        return DO_WOLFF
    else
        return DO_LOCAL
    end
end
