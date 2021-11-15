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
