struct CurieWeiss{T<:Integer, A<:AbstractVector{T}}
    spins::A
end

function CurieWeiss(N::Int, ::Type{T} = Int8) where {T<:Signed}
    spins = rand((one(T), -one(T)), N)
    return CurieWeiss(spins)
end

function CurieWeiss(N::Int, ::Type{Bool})
    spins = bitrand(N)
    return CurieWeiss(spins)
end

magnetization(model::CurieWeiss) = sum(ising, model.spins)

function energy(model::CurieWeiss, h::Real = false; f = nothing)
    N = length(model.spins)
    M = magnetization(model)
    E = 1/2 * (1 - M^2 / N) - h * M
    if isnothing(f)
        return E
    else
        return E + f(M)
    end
end

random_configuration!(model::CurieWeiss{<:Bool}) = rand!(model.spins)

function random_configuration!(model::CurieWeiss{<:Signed})
    return rand!(model.spins, (one(eltype(model.spins)), -one(eltype(model.spins))))
end

"""
    metropolis!(curie_wiess_model, β, h=0; f=nothing, steps=1, save_interval=length(spins))

Like `metropolis!`, but takes an additional argument `f`, which is a function of
the magnetization that gets added to the energy, `E = (...) + f(M)`.
"""
function metropolis!(model::CurieWeiss,
    β::Real,
    h::Real = false;
    f = nothing,
    steps::Int = 1,
    save_interval::Int = length(model.spins),
)
    @assert steps ≥ 1
    @assert save_interval ≥ 1

    #= Track history of magnetization and energy =#
    M0 = magnetization(model)
    if isnothing(f)
        E0 = energy(model, h)
    else
        E0 = energy(model, h) + f(M0)
    end

    M = zeros(typeof(M0), steps)
    E = zeros(typeof(E0), steps)

    M[1] = M0
    E[1] = E0

    #= Track the history of configurations only every 'save_interval' steps. =#
    spins_t = zeros(eltype(model.spins), length(model.spins), length(1:save_interval:steps))
    spins_t[:,1] .= model.spins

    for t ∈ 2:steps
        metropolis_step!(model, h; t = t, M = M, E = E, β = β, f = f)
        if t ∈ 1:save_interval:steps
            spins_t[:, cld(t, save_interval)] .= model.spins
        end
    end

    return spins_t, M, E
end

function metropolis_step!(model::CurieWeiss, h::Real = false; t, M, E, β, f = nothing)
    N = length(model.spins)
    i = rand(1:N)
    M[t] = M[t - 1]
    E[t] = E[t - 1]
    ΔM = -2 * ising(model.spins[i])
    ΔE = -(M[t] + ΔM/2) * ΔM / N + h * ΔM
    if !isnothing(f)
        Δf = f(M[t] + ΔM) - f(M[t])
        ΔE += Δf
    end
    if ΔE ≤ 0 || rand() < exp(-β * ΔE)
        model.spins[i] = flip(model.spins[i])
        M[t] += ΔM
        E[t] += ΔE
        return true
    else
        return false
    end
end
