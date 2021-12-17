"""
    metropolis_f!(spins, β, h = 0; f, steps = 1, save_interval = length(spins))

Like `metropolis!`, but takes an additional argument `f`, which is a function of
the magnetization that gets added to the energy, `E = (...) + f(M)`.
"""
function metropolis_f!(
    spins::AbstractMatrix,
    β::Real,
    h::Real = false;
    f,
    steps::Int = 1,
    save_interval::Int = length(spins),
)
    @assert steps ≥ 1
    @assert save_interval ≥ 1

    #= Track history of magnetization and energy =#
    M0 = magnetization(spins)
    E0 = energy(spins, h) + f(M0)

    M = zeros(typeof(M0), steps)
    E = zeros(typeof(E0), steps)

    M[1] = M0
    E[1] = E0

    #= Track the history of configurations only every 'save_interval' steps. =#
    spins_t = zeros(eltype(spins), size(spins)..., length(1:save_interval:steps))
    spins_t[:,:,1] .= spins

    for t ∈ 2:steps
        metropolis_step_f!(spins, h; t = t, M = M, E = E, β = β, f = f)
        if t ∈ 1:save_interval:steps
            spins_t[:, :, cld(t, save_interval)] .= spins
        end
    end

    return spins_t, M, E
end

function metropolis_step_f!(s::AbstractMatrix, h::Real = false; t, M, E, β, f)
    i, j = rand.(Base.OneTo.(size(s)))
    S = neighbor_sum_div_2(s, i, j)
    M[t] = M[t - 1]
    E[t] = E[t - 1]
    ΔM = -2 * ising(s[i,j])
    Δf = f(M[t - 1] + ΔM) - f(M[t - 1])
    ΔE = -ΔM * (2S + h) + Δf
    if ΔE ≤ 0 || randexp() > β * ΔE
        s[i,j] = flip(s[i,j])
        M[t] += ΔM
        E[t] += ΔE
        return true
    else
        return false
    end
end
