"""
    hybrid!(σ, β; steps = 1, save_interval = length(σ), local_steps = length(σ))

Hybrid sampler, performing a Wolff cluster move every `local_steps` steps, and Metropolis
steps otherwise. That is, `local_steps - 1` Metropolis steps, then one Wolff cluster move,
then another `local_steps - 1` Metropolis steps, one more Wolff cluster move, and so on.
"""
function hybrid!(
    σ::IsingArray,
    β::Real;
    steps::Int = 1,
    save_interval::Int = length(σ),
    local_steps::Int = length(σ) # one Wolff step every local_steps steps
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
    selectdim(σ_t, ndims(σ) + 1, 1) .= σ

    Padd = wolff_padd(β)
    queued = falses(size(σ)) # Wolff work buffer, reused across steps

    for t ∈ 2:steps
        if t ∈ 1:local_steps:steps # Wolff step
            wolff_step!(σ; t = t, M = M, E = E, Padd = Padd, queued = queued)
        else # Metropolis step
            metropolis_step!(σ; t = t, M = M, E = E, β = β)
        end
        if t ∈ 1:save_interval:steps
            selectdim(σ_t, ndims(σ) + 1, cld(t, save_interval)) .= σ
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
Each step is a Wolff move with probability proportional to the rate of spin flips per
second measured so far for Wolff moves (relative to Metropolis), and a Metropolis step
otherwise. Since this depends on timings, runs are not reproducible by seeding the RNG.
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
    selectdim(σ_t, ndims(σ) + 1, 1) .= σ

    Padd = wolff_padd(β)
    queued = falses(size(σ)) # Wolff work buffer, reused across steps

    for t ∈ 2:steps
        if hybrid_decide(hybrid_stats)
            hybrid_stats.wolff_time += @elapsed begin
                flipped = wolff_step!(σ; t = t, M = M, E = E, Padd = Padd, queued = queued)
                hybrid_stats.wolff_flip += flipped
            end
        else
            hybrid_stats.local_time += @elapsed begin
                flipped = metropolis_step!(σ; t = t, M = M, E = E, β = β)
                hybrid_stats.local_flip += flipped
            end
        end
        if t ∈ 1:save_interval:steps
            selectdim(σ_t, ndims(σ) + 1, cld(t, save_interval)) .= σ
        end
    end

    return σ_t, M, E
end

function hybrid_decide(hybrid_stats::HybridStats)
    # return true -> do Wolff, else do Metropolis
    DO_WOLFF = true
    DO_LOCAL = false

    # at the start, until both moves have been timed and have flipped some spin
    if iszero(hybrid_stats.wolff_time) || iszero(hybrid_stats.wolff_flip)
        return DO_WOLFF
    elseif iszero(hybrid_stats.local_time) || iszero(hybrid_stats.local_flip)
        return DO_LOCAL
    end

    #= The choice must not depend on the recent history of accepted moves, which depends on
    the current state. E.g., forcing Metropolis steps until one is accepted (to keep the
    flip counts balanced) biases the sampled distribution towards ordered states. =#
    wolff_rate = hybrid_stats.wolff_flip / hybrid_stats.wolff_time
    local_rate = hybrid_stats.local_flip / hybrid_stats.local_time

    Pwolff = wolff_rate / (wolff_rate + local_rate)
    if rand() < Pwolff
        return DO_WOLFF
    else
        return DO_LOCAL
    end
end
