"""
    neighbor_sum_matrix(N)

Returns a matrix `K`, such that (K * A)[i,j] = A[i+1,j] + A[i-1,j]
and (A * K)[i,j] = A[i,j+1] + A[i,j-1], up to periodic boundaries.
"""
neighbor_sum_matrix(N::Int) = I(N)[[2:N..., 1], :] + I(N)[[N, 1:(N - 1)...], :]

function adjacency_matrix(N::Int, M::Int)
    A = falses(N, M, N, M)
    for i in 1:N, j in 1:M, d in (-1, 1)
        A[i, j, mod1(i + d, N), j] = true
        A[i, j, i, mod1(j + d, M)] = true
    end
    return A
end

function energy(spins::AbstractMatrix)
    J = adjacency_matrix(size(spins)...)
    return vec(spins)' * J * vec(spins)
end

function metropolis!(
    spins::Matrix{Int8},
    β::Real;
    time::Int,
    save_interval = min(size(spins)...)
)
    #= We only need to evalute exp.(-β .* ΔE) for ΔE = 0, 1, 2, 3.
    Therefore we store a look-up table. =#
    _Expβ = exp.(-β .* 0:3)

    #= Track history of magnetization and energy =#
    magnetizations = zeros(time)
    energies = zeros(time)

    magnetizations[1] = mean(spins) # magnetization
    energies[1] = energy(spins) / length(spins) # energy per spin

    #= Track the history of configurations only every 'save_interval' steps. =#
    spins_t = [copy(spins)]

    for t = 2:time
        i, j = rand(1:N), rand(1:M)

        neighbors_sum = 0
        for d in (-1, 1)
            neighbors_sum += spins[mod1(i + d, N), j]
            neighbors_sum += spins[i, mod1(j + d, M)]
        end

        ΔE = -2 * spins[i,j] * neighbors_sum
        if ΔE < 0 || rand() < _Expβ[ΔE + 1]
            spins[i,j] = -spins[i,j]
            magnetizations[t] = magnetizations[t - 1] + 2spins[i,j] / length(spins)
            energies[t] = energies[t - 1] + ΔE
        end

        if save_every !== nothing && t % save_every == 0
            push!(spins_t, copy(spins))
        end
    end

    return spins_t, magnetizations, energies
end
