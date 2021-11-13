# Examples

## [Example simulation with Metropolis](@id example_metropolis)

```@example
using SquareIsingModel: metropolis!, onsager_magnetization, βc
using Statistics, CairoMakie, Random, ProgressMeter

Random.seed!(1) # make reproducible

βs = 0:0.05:1 # inverse temperatures to simulate
mavg = zeros(length(βs))
mstd = zeros(length(βs))
spins = rand((Int8(1), Int8(-1)), 50, 50);
@showprogress for (k,β) in enumerate(βs)
    spins_t, m, E = metropolis!(spins, β, 10^7)
    mavg[k] = mean(m[(length(m) ÷ 2):end])
    mstd[k] = std(m[(length(m) ÷ 2):end])
end
mavg .= abs.(mavg) # remove sign invariance

fig = Figure(resolution=(600,400))
ax = Axis(fig[1,1], xlabel="β", ylabel="m")
lines!(ax, 0:0.01:1, onsager_magnetization, color=:black, label="analytical")
scatter!(ax, βs, mavg, color=:red, markersize=5, label="MC")
errorbars!(ax, βs, mavg, mstd/2, color=:red, whiskerwidth=5)
axislegend(ax, position=:rb)
fig
```

## [Example simulation with Wolff](@id example_wolff)

```@example
using SquareIsingModel: metropolis!, onsager_magnetization, βc
using Statistics, CairoMakie, Random, ProgressMeter

Random.seed!(1) # make reproducible

βs = 0:0.05:1 # inverse temperatures to simulate
mavg = zeros(length(βs))
mstd = zeros(length(βs))
spins = rand((Int8(1), Int8(-1)), 50, 50);
@showprogress for (k,β) in enumerate(βs)
    spins_t, m, E = wolff!(spins, β, 10^3)
    mavg[k] = mean(m[(length(m) ÷ 2):end])
    mstd[k] = std(m[(length(m) ÷ 2):end])
end
mavg .= abs.(mavg) # remove sign invariance

fig = Figure(resolution=(600,400))
ax = Axis(fig[1,1], xlabel="β", ylabel="m")
lines!(ax, 0:0.01:1, onsager_magnetization, color=:black, label="analytical")
scatter!(ax, βs, mavg, color=:red, markersize=5, label="MC")
errorbars!(ax, βs, mavg, mstd/2, color=:red, whiskerwidth=5)
axislegend(ax, position=:rb)
fig
```
