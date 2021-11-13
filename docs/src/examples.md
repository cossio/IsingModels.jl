# Examples

## [Example simulation with Metropolis](@id example_metropolis)

```@example
import SquareIsingModel as Ising
using Statistics, CairoMakie, Random, ProgressMeter

Random.seed!(1) # make reproducible

βs = 0:0.05:1 # inverse temperatures to simulate
mavg = zeros(length(βs))
mstd = zeros(length(βs))
spins = Ising.random_configuration(50, 50)
@showprogress for (k,β) in enumerate(βs)
    spins_t, m, E = Ising.metropolis!(spins, β, 10^7)
    mavg[k] = mean(m[(length(m) ÷ 2):end])
    mstd[k] = std(m[(length(m) ÷ 2):end])
end
mavg .= abs.(mavg) # remove sign invariance

fig = Figure(resolution=(600,400))
ax = Axis(fig[1,1], xlabel="β", ylabel="m")
lines!(ax, 0:0.01:1, Ising.onsager_magnetization, color=:black, label="analytical")
scatter!(ax, βs, mavg, color=:red, markersize=5, label="MC")
errorbars!(ax, βs, mavg, mstd/2, color=:red, whiskerwidth=5)
axislegend(ax, position=:rb)
fig
```

## [Example simulation with Wolff](@id example_wolff)

```@example
import SquareIsingModel as Ising
using Statistics, CairoMakie, Random, ProgressMeter

Random.seed!(1) # make reproducible

βs = 0:0.05:1 # inverse temperatures to simulate
mavg = zeros(length(βs))
mstd = zeros(length(βs))
spins = Ising.random_configuration(50, 50)
@showprogress for (k,β) in enumerate(βs)
    spins_t, m, E = Ising.wolff!(spins, β, 10^3)
    mavg[k] = mean(m[(length(m) ÷ 2):end])
    mstd[k] = std(m[(length(m) ÷ 2):end])
end
mavg .= abs.(mavg) # remove sign invariance

fig = Figure(resolution=(600,400))
ax = Axis(fig[1,1], xlabel="β", ylabel="m")
lines!(ax, 0:0.01:1, Ising.onsager_magnetization, color=:black, label="analytical")
scatter!(ax, βs, mavg, color=:red, markersize=5, label="MC")
errorbars!(ax, βs, mavg, mstd/2, color=:red, whiskerwidth=5)
axislegend(ax, position=:rb)
fig
```

## Example Wolff clusters

```@example
import SquareIsingModel as Ising
using Random, Colors, ColorSchemes, CairoMakie

Random.seed!(69) # reproducibility

β = 0.6
spins = Ising.random_configuration(50, 50)
Ising.metropolis!(spins, β, 10^7)
cluster = Ising.wolff_cluster(spins, 25, 25, 1 - exp(-2β))

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