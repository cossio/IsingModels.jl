"""
    metropolis!(σ, β, h = 0; steps = 1, save_interval = length(σ))

Performs one or more Metropolis MC steps from the configuration `σ`, at inverse
temperature `β`.
Returns three lists: `σ_t, M, E`, where `σ_t` contains configurations
sampled at intervals `save_interval` (by default equals the number of sites),
`M` is the record of magnetizations, and `E` the record of energies.
"""
function metropolis!(
    σ::IsingArray,
    β::Real,
    h::Real = false;
    steps::Int = 1,
    save_interval::Int = length(σ),
    f = nothing
)
    @assert steps ≥ 1
    @assert save_interval ≥ 1

    #= Track history of magnetization and energy =#
    M0 = magnetization(σ)
    M = zeros(typeof(M0), steps)
    M[1] = M0

    E0 = energy(σ, h; f = f)
    E = zeros(typeof(E0), steps)
    E[1] = E0

    #= Track the history of configurations only every 'save_interval' steps. =#
    σ_t = falses(size(σ)..., length(1:save_interval:steps))
    selectdim(σ_t, ndims(σ) + 1, 1) .= σ

    for t ∈ 2:steps
        metropolis_step!(σ, h; t = t, M = M, E = E, β = β, f = f)
        if t ∈ 1:save_interval:steps
            selectdim(σ_t, ndims(σ) + 1, cld(t, save_interval)) .= σ
        end
    end

    return σ_t, M, E
end

function metropolis_step!(σ::IsingArray, h::Real = false; t, M, E, β, f = nothing)
    i = rand(CartesianIndices(σ))
    S = sum_neighbors(σ, i)
    M[t] = M[t - 1]
    E[t] = E[t - 1]
    ΔM = -2spin(σ[i])

    if isnothing(f)
        Δf = 0
    else
        Δf = f(M[t] + ΔM) - f(M[t])
    end

    ΔE = -ΔM * (S + h) + Δf

    if ΔE ≤ 0 || randexp() > β * ΔE
        σ[i] = flip(σ[i])
        M[t] += ΔM
        E[t] += ΔE
        return true
    else
        return false
    end
end
