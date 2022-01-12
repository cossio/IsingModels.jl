#=
# Examples using Wolff cluster sampling method

## Magnetization as a function of temperature

First load some packages
=#

import IsingModels as Ising
using Statistics, CairoMakie, Random

Random.seed!(1) # make reproducible
nothing #hide

#=
First example.
=#

βs = 0:0.025:1
fig = Figure(resolution=(600,400))
ax = Axis(fig[1,1], xlabel="β", ylabel="m")
lines!(ax, 0:0.01:1, Ising.onsager_magnetization, color=:black, label="analytical")

mavg = zeros(length(βs))
mstd = zeros(length(βs))
σ = bitrand(32, 32)
for (k,β) in enumerate(βs)
    σ_t, M, E = Ising.wolff!(σ, β, steps=10^3)
    m = abs.(M[(length(M) ÷ 2):end]) / length(σ)
    mavg[k] = mean(m)
    mstd[k] = std(m)
end
scatter!(ax, βs, mavg, color=:blue, markersize=5, label="MC, L=30")
errorbars!(ax, βs, mavg, mstd/2, color=:blue, whiskerwidth=5)

mavg = zeros(length(βs))
mstd = zeros(length(βs))
σ = bitrand(64, 64)
for (k,β) in enumerate(βs)
    σ_t, M, E = Ising.wolff!(σ, β, steps=10^3)
    m = abs.(M[(length(M) ÷ 2):end]) / length(σ)
    mavg[k] = mean(m)
    mstd[k] = std(m)
end
scatter!(ax, βs, mavg, color=:red, markersize=5, label="MC, L=70")
errorbars!(ax, βs, mavg, mstd/2, color=:red, whiskerwidth=5)

axislegend(ax, position=:rb)
fig

#=
## Typical Wolff clusters at criticality
=#

import IsingModels as Ising
using Random, Colors, ColorSchemes, CairoMakie

Random.seed!(72) # reproducibility

β = Ising.βc
σ = bitrand(512, 512)
Ising.metropolis!(σ, β; steps=10^7)
Ising.wolff!(σ, β, steps=200)

cluster = Ising.wolff_cluster(σ, CartesianIndex(256, 256), Ising.wolff_padd(β))

fig = Figure(resolution=(700, 300))

ax = Axis(fig[1,1], title="spins")
hmap = heatmap!(ax, σ, colormap=cgrad([:purple, :orange], [0.5]; categorical=true))
cbar = Colorbar(fig[1,2], hmap)
cbar.ticks = ([-0.5, 0.5], ["-1", "1"])

ax = Axis(fig[1,3], title="cluster")
hmap = heatmap!(ax, cluster, colormap=cgrad([:white, :black], [0.5]; categorical=true))
cbar = Colorbar(fig[1, 4], hmap)
cbar.ticks = ([0.25, 0.75], ["0", "1"])

fig

#=
## Average size of Wolff's clusters as a function of temperature
=#

import IsingModels as Ising
using Random, Statistics, Colors, ColorSchemes, CairoMakie

Random.seed!(1) # make reproducible

Ts = 0:0.2:5
βs = inv.(Ts)
fig = Figure(resolution=(600,400))
ax = Axis(fig[1,1], xlabel=L"temperature $T$", ylabel=L"average Wolff's cluster size / $N$")

clavg = zeros(length(βs))
clstd = zeros(length(βs))
σ = bitrand(32, 32)
for (k,β) in enumerate(βs)
    σ_t, M, E = Ising.wolff!(σ, β, steps=10^3)
    cluster_size = abs.(M[2:end] - M[1:(end - 1)]) .÷ 2
    clavg[k] = mean(cluster_size / length(σ))
    clstd[k] = std(cluster_size / length(σ))
end
scatter!(ax, Ts, clavg, color=:blue, markersize=5, label="L=30")
lines!(ax, Ts, clavg, color=:blue)
errorbars!(ax, Ts, clavg, clstd/2, color=:blue, whiskerwidth=5)

clavg = zeros(length(βs))
clstd = zeros(length(βs))
σ = bitrand(64, 64)
for (k,β) in enumerate(βs)
    σ_t, M, E = Ising.wolff!(σ, β, steps=10^3)
    cluster_size = abs.(M[2:end] - M[1:(end - 1)]) .÷ 2
    clavg[k] = mean(cluster_size / length(σ))
    clstd[k] = std(cluster_size / length(σ))
end
scatter!(ax, Ts, clavg, color=:red, markersize=5, label="L=70")
lines!(ax, Ts, clavg, color=:red)
errorbars!(ax, Ts, clavg, clstd/2, color=:red, whiskerwidth=5)

axislegend(ax, position=:rt)
fig

#=
## Wolff explores configurations efficiently at the critical temperature

The following example shows how at the critical temperature, the Metropolis sampler gets stuck in particular cluster structures.
In contrast the Wolff algorihm explores diverse states.
=#

using Statistics, CairoMakie, Random
import IsingModels as Ising

Random.seed!(3) # reproducibility

L = 128
N = L^2
σ_t_metro, M_metro, E_metro = Ising.metropolis!(falses(L, L), Ising.βc; steps=10^6, save_interval=2*10^5)
σ_t_wolff, M_wolff, E_wolff = Ising.wolff!(falses(L, L), Ising.βc, steps=10^4, save_interval=2*10^3)

fig = Figure(resolution=(1000, 650))
for t in 1:size(σ_t_metro, 3)
    ax = Axis(fig[1,t], title="t=$t, metropolis")
    heatmap!(ax, σ_t_metro[:,:,t], colormap=cgrad([:purple, :orange], [0.5]; categorical=true))
    hidedecorations!(ax)
end
for t in 1:size(σ_t_wolff, 3)
    ax = Axis(fig[2,t], title="t=$t, wolff")
    heatmap!(ax, σ_t_wolff[:,:,t], colormap=cgrad([:purple, :orange], [0.5]; categorical=true))
    hidedecorations!(ax)
end

ax = Axis(fig[3,1:2], title="magnetization")
lines!(ax, M_wolff[1:10:end] / N, label="wolff")
lines!(ax, M_metro[1:1000:end] / N, label="metropolis", linewidth=2, color=:red)
ylims!(ax, (-1,1))
ax = Axis(fig[3,3:4], title="energy")
lines!(ax, E_wolff[1:10:end] / N, label="wolff")
lines!(ax, E_metro[1:1000:end] / N, label="metropolis", linewidth=2, color=:red)
axislegend(ax, position = :rb)

fig


#=
## Binder's parameter

As the system size grows, the crossing point of the different curves is the critical temperature.
=#

using Statistics, CairoMakie, Random
import IsingModels as Ising

Random.seed!(1) # make reproducible
Ts = 2.2:0.01:2.3
βs = inv.(Ts)

fig = Figure(resolution=(600, 400))
ax = Axis(fig[1,1], xlabel=L"temperature $T$ ($=1/\beta$)", ylabel="Binder parameter")
@time for (L, color) in zip([4, 8, 16, 32], [:green, :orange, :blue, :red])
    U = zeros(length(βs))
    σ = bitrand(L, L)
    for (k, β) in enumerate(βs)
        σ_t, M, E = Ising.wolff!(σ, β, steps=10^5)
        U[k] = (3 - mean(M.^4) / mean(M.^2)^2) / 2
    end
    scatter!(ax, Ts, U, color=color, markersize=5, label=L"L=%$L")
    lines!(ax, Ts, U, color=color, markersize=5)
end
vlines!(ax, [1 / Ising.βc], label=L"Onsager's $T_c$", color=:black, linestyle=:dash)
axislegend(ax, position=:lb)

fig

#=
## Heat capacity vs. exact expression
=#

using Statistics, CairoMakie, Random
import IsingModels as Ising

Random.seed!(1) # make reproducible
Ts = 1.8:0.01:3
βs = inv.(Ts)

fig = Figure(resolution=(800, 400))
ax = Axis(fig[1,1], xlabel=L"temperature $T$ ($=1/\beta$)", ylabel="heat capacity", limits=(extrema(Ts)..., 0,2))

@time for (L, color) in zip([4, 8, 16, 32], [:green, :orange, :blue, :red])
    C = zeros(length(βs))
    for (k, β) in enumerate(βs)
        σ = bitrand(L, L)
        σ_t, M, E = Ising.wolff!(σ, β, steps=10^5)
        C[k] = β^2/length(σ) * var(E)
    end
    scatter!(ax, Ts, C, color=color, markersize=5, label=L"L=%$L")
end
lines!(ax, Ts, Ising.onsager_heat_capacity.(βs), color=:black, label="exact")
vlines!(ax, [1 / Ising.βc], label=L"Onsager's $T_c$", color=:black, linestyle=:dash)
Legend(fig[1,2], ax)

fig

#=
## Internal energy vs. exact expression
=#

using Statistics, CairoMakie, Random
import IsingModels as Ising

Random.seed!(1) # make reproducible
Ts = 1.8:0.05:3
βs = inv.(Ts)

fig = Figure(resolution=(800, 400))
ax = Axis(fig[1,1], xlabel=L"temperature $T$ ($=1/\beta$)", ylabel="internal energy")

@time for (L, color) in zip([4, 8, 16, 32], [:green, :orange, :blue, :red])
    Eavg = zeros(length(βs))
    Estd = zeros(length(βs))
    for (k, β) in enumerate(βs)
        σ = bitrand(L, L)
        σ_t, M, E = Ising.wolff!(σ, β, steps=10^5)
        Eavg[k] = mean(E / length(σ))
        Estd[k] = std(E / length(σ))
    end
    scatter!(ax, Ts, Eavg, color=color, markersize=5, label=L"L=%$L")
    errorbars!(ax, Ts, Eavg, Estd/2, color=color, whiskerwidth=5)
end
lines!(ax, Ts, Ising.onsager_internal_energy.(βs), color=:black, label="exact")
vlines!(ax, [1 / Ising.βc], label=L"Onsager's $T_c$", color=:black, linestyle=:dash)
Legend(fig[1,2], ax)

fig
