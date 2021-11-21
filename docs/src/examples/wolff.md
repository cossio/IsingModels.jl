# Examples using Wolff cluster sampling method


## Magnetization as a function of temperature

```@example
import SquareIsingModel as Ising
using Statistics, CairoMakie, Random, ProgressMeter

Random.seed!(1) # make reproducible

βs = 0:0.025:1
fig = Figure(resolution=(600,400))
ax = Axis(fig[1,1], xlabel="β", ylabel="m")
lines!(ax, 0:0.01:1, Ising.onsager_magnetization, color=:black, label="analytical")

mavg = zeros(length(βs))
mstd = zeros(length(βs))
spins = Ising.random_configuration(30)
@showprogress for (k,β) in enumerate(βs)
    spins_t, M, E = Ising.wolff!(spins, β, 10^3)
    m = abs.(M[(length(M) ÷ 2):end]) / length(spins)
    mavg[k] = mean(m)
    mstd[k] = std(m)
end
scatter!(ax, βs, mavg, color=:blue, markersize=5, label="MC, L=30")
errorbars!(ax, βs, mavg, mstd/2, color=:blue, whiskerwidth=5)

mavg = zeros(length(βs))
mstd = zeros(length(βs))
spins = Ising.random_configuration(70)
@showprogress for (k,β) in enumerate(βs)
    spins_t, M, E = Ising.wolff!(spins, β, 10^3)
    m = abs.(M[(length(M) ÷ 2):end]) / length(spins)
    mavg[k] = mean(m)
    mstd[k] = std(m)
end
scatter!(ax, βs, mavg, color=:red, markersize=5, label="MC, L=70")
errorbars!(ax, βs, mavg, mstd/2, color=:red, whiskerwidth=5)

axislegend(ax, position=:rb)
fig
```


## Typical Wolff clusters at criticality

```@example
import SquareIsingModel as Ising
using Random, Colors, ColorSchemes, CairoMakie

Random.seed!(62) # reproducibility

β = Ising.βc
spins = Ising.random_configuration(400)
Ising.metropolis!(spins, β, 10^7)
Ising.wolff!(spins, β, 200)

cluster = Ising.wolff_cluster(spins, 200, 200, Ising.wolff_padd(β))

fig = Figure(resolution=(700, 300))

ax = Axis(fig[1,1], title="spins")
hmap = heatmap!(ax, spins, colormap=cgrad([:purple, :orange], [0.5]; categorical=true))
cbar = Colorbar(fig[1,2], hmap)
cbar.ticks = ([-0.5, 0.5], ["-1", "1"])

ax = Axis(fig[1,3], title="cluster")
hmap = heatmap!(ax, cluster, colormap=cgrad([:white, :black], [0.5]; categorical=true))
cbar = Colorbar(fig[1, 4], hmap)
cbar.ticks = ([0.25, 0.75], ["0", "1"])

fig
```


## Average size of Wolff's clusters as a function of temperature

```@example
import SquareIsingModel as Ising
using Random, Statistics, Colors, ColorSchemes, CairoMakie, ProgressMeter

Random.seed!(1) # make reproducible

Ts = 0:0.2:5
βs = inv.(Ts)
fig = Figure(resolution=(600,400))
ax = Axis(fig[1,1], xlabel=L"temperature $T$", ylabel=L"average Wolff's cluster size / $N$")

clavg = zeros(length(βs))
clstd = zeros(length(βs))
spins = Ising.random_configuration(30)
@showprogress for (k,β) in enumerate(βs)
    spins_t, M, E = Ising.wolff!(spins, β, 10^3)
    cluster_size = abs.(M[2:end] - M[1:(end - 1)]) .÷ 2
    clavg[k] = mean(cluster_size / length(spins))
    clstd[k] = std(cluster_size / length(spins))
end
scatter!(ax, Ts, clavg, color=:blue, markersize=5, label="L=30")
lines!(ax, Ts, clavg, color=:blue)
errorbars!(ax, Ts, clavg, clstd/2, color=:blue, whiskerwidth=5)

clavg = zeros(length(βs))
clstd = zeros(length(βs))
spins = Ising.random_configuration(70)
@showprogress for (k,β) in enumerate(βs)
    spins_t, M, E = Ising.wolff!(spins, β, 10^3)
    cluster_size = abs.(M[2:end] - M[1:(end - 1)]) .÷ 2
    clavg[k] = mean(cluster_size / length(spins))
    clstd[k] = std(cluster_size / length(spins))
end
scatter!(ax, Ts, clavg, color=:red, markersize=5, label="L=70")
lines!(ax, Ts, clavg, color=:red)
errorbars!(ax, Ts, clavg, clstd/2, color=:red, whiskerwidth=5)

axislegend(ax, position=:rt)
fig
```


## Wolff explores configurations efficiently at the critical temperature

The following example shows how at the critical temperature, the Metropolis sampler gets stuck in particular cluster structures.
In contrast the Wolff algorihm explores diverse states.


```@example
using Statistics, CairoMakie, Random
import SquareIsingModel as Ising

Random.seed!(3) # reproducibility

L = 100
N = L^2
spins = Ising.random_configuration(100)
spins .= 1
spins_t_metro, M_metro, E_metro = Ising.metropolis!(spins, Ising.βc, 10^6);
spins .= 1
spins_t_wolff, M_wolff, E_wolff = Ising.wolff!(spins, Ising.βc, 10^4);

fig = Figure(resolution=(1000, 650))
for (col, t) in enumerate(1:20:100)
    ax = Axis(fig[1,col], title="t=$t, metropolis")
    heatmap!(ax, spins_t_metro[:,:,t], colormap=cgrad([:purple, :orange], [0.5]; categorical=true))
    hidedecorations!(ax)
end
for (col, t) in enumerate(1:2000:10000)
    ax = Axis(fig[2,col], title="t=$t, wolff")
    heatmap!(ax, spins_t_wolff[:,:,t], colormap=cgrad([:purple, :orange], [0.5]; categorical=true))
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
```


## Binder's parameter

As the system size grows, the crossing point of the different curves is the critical temperature.

```@example
using Statistics, CairoMakie, Random
import SquareIsingModel as Ising

Random.seed!(1) # make reproducible
Ts = 2.2:0.01:2.3
βs = inv.(Ts)

fig = Figure(resolution=(600, 400))
ax = Axis(fig[1,1], xlabel=L"temperature $T$ ($=1/\beta$)", ylabel="Binder parameter")
@time for (L, color) in zip([4, 8, 16, 32], [:green, :orange, :blue, :red])
    U = zeros(length(βs))
    spins = Ising.random_configuration(L)
    for (k, β) in enumerate(βs)
        spins_t, M, E = Ising.wolff!(spins, β, 10^5)
        U[k] = (3 - mean(M.^4) / mean(M.^2)^2) / 2
    end
    scatter!(ax, Ts, U, color=color, markersize=5, label=L"L=%$L")
    lines!(ax, Ts, U, color=color, markersize=5)
end
vlines!(ax, [1 / Ising.βc], label=L"Onsager's $T_c$", color=:black, linestyle=:dash)
axislegend(ax, position=:lb) 

fig
```


## Heat capacity vs. exact expression

```@example
using Statistics, CairoMakie, Random
import SquareIsingModel as Ising

Random.seed!(1) # make reproducible
Ts = 1.8:0.01:3
βs = inv.(Ts)

fig = Figure(resolution=(800, 400))
ax = Axis(fig[1,1], xlabel=L"temperature $T$ ($=1/\beta$)", ylabel="heat capacity", limits=(extrema(Ts)..., 0,2))

@time for (L, color) in zip([4, 8, 16, 32], [:green, :orange, :blue, :red])
    C = zeros(length(βs))
    for (k, β) in enumerate(βs)
        spins = Ising.random_configuration(L)
        spins_t, M, E = Ising.wolff!(spins, β, 10^5)
        C[k] = β^2/length(spins) * var(E)
    end
    scatter!(ax, Ts, C, color=color, markersize=5, label=L"L=%$L")
end
lines!(ax, Ts, Ising.onsager_heat_capacity.(βs), color=:black, label="exact")
vlines!(ax, [1 / Ising.βc], label=L"Onsager's $T_c$", color=:black, linestyle=:dash)
Legend(fig[1,2], ax)

fig
```


## Internal energy vs. exact expression

```@example
using Statistics, CairoMakie, Random
import SquareIsingModel as Ising

Random.seed!(1) # make reproducible
Ts = 1.8:0.05:3
βs = inv.(Ts)

fig = Figure(resolution=(800, 400))
ax = Axis(fig[1,1], xlabel=L"temperature $T$ ($=1/\beta$)", ylabel="internal energy")

@time for (L, color) in zip([4, 8, 16, 32], [:green, :orange, :blue, :red])
    Eavg = zeros(length(βs))
    Estd = zeros(length(βs))
    for (k, β) in enumerate(βs)
        spins = Ising.random_configuration(L)
        spins_t, M, E = Ising.wolff!(spins, β, 10^5)
        Eavg[k] = mean(E / length(spins))
        Estd[k] = std(E / length(spins))
    end
    scatter!(ax, Ts, Eavg, color=color, markersize=5, label=L"L=%$L")
    errorbars!(ax, Ts, Eavg, Estd/2, color=color, whiskerwidth=5)
end
lines!(ax, Ts, Ising.onsager_internal_energy.(βs), color=:black, label="exact")
vlines!(ax, [1 / Ising.βc], label=L"Onsager's $T_c$", color=:black, linestyle=:dash)
Legend(fig[1,2], ax)

fig
```
