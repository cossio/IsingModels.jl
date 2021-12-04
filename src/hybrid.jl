"""
    hybrid!(spins, β; steps = 1, save_interval = length(spins), local_steps = length(spins))

Hybrid sampler, performing `local_steps` of Metropolis sampling, then one Wolff cluster
move, then another `local_steps` of Metropolis sampling, one more Wolff cluster move,
and so on.
"""
function hybrid!(
    spins::AbstractMatrix,
    β::Real;
    steps::Int = 1,
    save_interval::Int = length(spins),
    local_steps::Int = length(spins) # Metropolis step per Wolff step
)
    @assert steps ≥ 1
    @assert save_interval ≥ 1
    @assert local_steps ≥ 1

    #= Track history of magnetization and energy =#
    M = zeros(Int, steps)
    E = zeros(Int, steps)

    M[1] = magnetization(spins)
    E[1] = energy(spins)

    #= Track the history of configurations only every 'save_interval' steps. =#
    spins_t = zeros(eltype(spins), size(spins)..., length(1:save_interval:steps))
    spins_t[:,:,1] .= spins

    Pmet = metropolis_acceptance_probabilities(β)
    Padd = wolff_padd(β)

    for t ∈ 2:steps
        if t ∈ 1:local_steps:steps # Wolff step
            wolff_step!(spins; t = t, M = M, E = E, Padd = Padd)
        else # Metropolis step
            metropolis_step!(spins; t = t, M = M, E = E, Paccept = Pmet)
        end
        if t ∈ 1:save_interval:steps
            spins_t[:, :, cld(t, save_interval)] .= spins
        end
    end

    return spins_t, M, E
end

mutable struct HybridStats
    local_flip::Float64
    wolff_flip::Float64
    local_time::Float64
    wolff_time::Float64
end
HybridStats() = HybridStats(0, 0, 0, 0)

"""
    dynamic_hybrid!(spins, β; steps, save_interval)

Same as `hybrid!`, but adjusts numbers of Metropolis and Wolff steps dynamically.
"""
function dynamic_hybrid!(
    spins::AbstractMatrix,
    β::Real;
    steps::Int = 1,
    save_interval::Int = length(spins),
    hybrid_stats::HybridStats = HybridStats()
)
    @assert steps ≥ 1
    @assert save_interval ≥ 1

    #= Track history of magnetization and energy =#
    M = zeros(Int, steps)
    E = zeros(Int, steps)

    M[1] = sum(spins) # magnetization
    E[1] = energy(spins)

    #= Track the history of configurations only every 'save_interval' steps. =#
    spins_t = zeros(eltype(spins), size(spins)..., length(1:save_interval:steps))
    spins_t[:,:,1] .= spins

    Pmet = metropolis_acceptance_probabilities(β)
    Padd = wolff_padd(β)

    for t ∈ 2:steps
        if hybrid_decide(hybrid_stats, length(spins))
            hybrid_stats.wolff_time += @elapsed begin
                flipped = wolff_step!(spins; t = t, M = M, E = E, Padd = Padd)
                hybrid_stats.wolff_flip += flipped
            end
        else
            hybrid_stats.local_time += @elapsed begin
                flipped = metropolis_step!(spins; t = t, M = M, E = E, Paccept = Pmet)
                hybrid_stats.local_flip += flipped
            end
        end
        if t ∈ 1:save_interval:steps
            spins_t[:, :, cld(t, save_interval)] .= spins
        end
    end

    return spins_t, M, E
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
