# Examples


## [Example simulation with Metropolis](@id example_metropolis)

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
spins = Ising.random_configuration(20)
@showprogress for (k,β) in enumerate(βs)
    spins_t, m, E = Ising.metropolis!(spins, β, 10^7)
    m_ = abs.(m[(length(m) ÷ 2):end])
    mavg[k] = mean(m_)
    mstd[k] = std(m_)
end
scatter!(ax, βs, mavg, color=:blue, markersize=5, label="MC, L=20")
errorbars!(ax, βs, mavg, mstd/2, color=:blue, whiskerwidth=5)

mavg = zeros(length(βs))
mstd = zeros(length(βs))
spins = Ising.random_configuration(50)
@showprogress for (k,β) in enumerate(βs)
    spins_t, m, E = Ising.metropolis!(spins, β, 10^7)
    m_ = abs.(m[(length(m) ÷ 2):end])
    mavg[k] = mean(m_)
    mstd[k] = std(m_)
end
scatter!(ax, βs, mavg, color=:red, markersize=5, label="MC, L=50")
errorbars!(ax, βs, mavg, mstd/2, color=:red, whiskerwidth=5)

axislegend(ax, position=:rb)
fig
```


## [Example simulation with Wolff](@id example_wolff)

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
    spins_t, m, E = Ising.wolff!(spins, β, 10^3)
    m_ = abs.(m[(length(m) ÷ 2):end])
    mavg[k] = mean(m_)
    mstd[k] = std(m_)
end
scatter!(ax, βs, mavg, color=:blue, markersize=5, label="MC, L=30")
errorbars!(ax, βs, mavg, mstd/2, color=:blue, whiskerwidth=5)

mavg = zeros(length(βs))
mstd = zeros(length(βs))
spins = Ising.random_configuration(70)
@showprogress for (k,β) in enumerate(βs)
    spins_t, m, E = Ising.wolff!(spins, β, 10^3)
    m_ = abs.(m[(length(m) ÷ 2):end])
    mavg[k] = mean(m_)
    mstd[k] = std(m_)
end
scatter!(ax, βs, mavg, color=:red, markersize=5, label="MC, L=70")
errorbars!(ax, βs, mavg, mstd/2, color=:red, whiskerwidth=5)

axislegend(ax, position=:rb)
fig
```


## Example Wolff clusters

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


## Average size of Wolff's clusters

```@example
import SquareIsingModel as Ising
using Random, Colors, ColorSchemes, CairoMakie

Random.seed!(1) # make reproducible

Ts = 0:0.2:5
βs = inv.(Ts)
fig = Figure(resolution=(600,400))
ax = Axis(fig[1,1], xlabel=L"temperature $T$", ylabel=L"average Wolff's cluster size / $N$")

clavg = zeros(length(βs))
clstd = zeros(length(βs))
spins = Ising.random_configuration(30)
@showprogress for (k,β) in enumerate(βs)
    spins_t, m, E = Ising.wolff!(spins, β, 10^3)
    cluster_size = abs.(m[2:end] - m[1:(end - 1)]) .* length(spins) / 2
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
    spins_t, m, E = Ising.wolff!(spins, β, 10^3)
    cluster_size = abs.(m[2:end] - m[1:(end - 1)]) .* length(spins) / 2
    clavg[k] = mean(cluster_size / length(spins))
    clstd[k] = std(cluster_size / length(spins))
end
scatter!(ax, Ts, clavg, color=:red, markersize=5, label="L=70")
lines!(ax, Ts, clavg, color=:red)
errorbars!(ax, Ts, clavg, clstd/2, color=:red, whiskerwidth=5)

axislegend(ax, position=:rt)
fig
```