"""
    hybrid!(spins, β, steps = 1; save_interval = length(spins), local_steps = length(spins))

Hybrid sampler, performing `local_steps` of Metropolis sampling, then one Wolff cluster
move, then another `local_steps` of Metropolis sampling, one more Wolff cluster move,
and so on.
"""
function hybrid!(
    spins::AbstractMatrix{Int8},
    β::Real, steps::Int = 1;
    save_interval::Int = length(spins),
    local_steps::Int = length(spins) # Metropolis step per Wolff step
)
    @assert steps ≥ 1
    @assert save_interval ≥ 1
    @assert local_steps ≥ 1

    #= Track history of magnetization and energy =#
    M = zeros(Int, steps)
    E = zeros(Int, steps)

    M[1] = sum(spins) # magnetization
    E[1] = energy(spins)

    #= Track the history of configurations only every 'save_interval' steps. =#
    spins_t = zeros(Int8, size(spins)..., length(1:save_interval:steps))
    spins_t[:,:,1] .= spins

    _Exp2β = metropolis_acceptance_probabilities(β)
    Padd = wolff_padd(β)

    for t ∈ 2:steps
        if t ∈ 1:local_steps:steps # Wolff step
            wolff_step!(spins; t = t, M = M, E = E, Padd = Padd)
        else # Metropolis step
            metropolis_step!(spins; t = t, M = M, E = E, _Exp2β = _Exp2β)
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
    dynamic_hybrid!(spins, β, steps; save_interval)

Same as `hybrid!`, but adjusts numbers of Metropolis and Wolff steps dynamically.
"""
function dynamic_hybrid!(
    spins::AbstractMatrix{Int8},
    β::Real, steps::Int = 1;
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
    spins_t = zeros(Int8, size(spins)..., length(1:save_interval:steps))
    spins_t[:,:,1] .= spins

    _Exp2β = metropolis_acceptance_probabilities(β)
    Padd = wolff_padd(β)

    wolff_flip = local_flip = 0.0
    wolff_time = local_time = 0.0

    for t ∈ 2:steps
        if hybrid_decide(wolff_flip, local_flip, wolff_time, local_time, length(spins))
            wolff_time += @elapsed begin
                wolff_flip += wolff_step!(spins; t = t, M = M, E = E, Padd = Padd)
            end
        else
            local_time += @elapsed begin
                local_flip += metropolis_step!(spins; t = t, M = M, E = E, _Exp2β = _Exp2β)
            end
        end
        if t ∈ 1:save_interval:steps
            spins_t[:, :, cld(t, save_interval)] .= spins
        end
    end

    hybrid_stats.local_flip = local_flip
    hybrid_stats.wolff_flip = wolff_flip
    hybrid_stats.local_time = local_time
    hybrid_stats.wolff_time = wolff_time

    return spins_t, M, E
end

function hybrid_decide(
    wolff_flip::Real,
    local_flip::Real,
    wolff_time::Real,
    local_time::Real,
    N::Integer
)
    # return true -> do Wolff, else do Metropolis
    DO_WOLFF = true
    DO_LOCAL = false

    if iszero(wolff_time) || iszero(wolff_flip) || wolff_flip * N < local_flip
        return DO_WOLFF
    elseif iszero(local_time) || iszero(local_flip) || local_flip * N < wolff_flip
        return DO_LOCAL
    end

    wolff_rate = wolff_flip / wolff_time
    local_rate = local_flip / local_time

    Pwolff = wolff_rate / (wolff_rate + local_rate)
    if rand() < Pwolff
        return DO_WOLFF
    else
        return DO_LOCAL
    end
end
