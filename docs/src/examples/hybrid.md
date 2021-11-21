# Examples using hybrid (Metropolis + Wolff) sampling method


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
    spins_t, M, E = Ising.hybrid!(spins, β, 10^7)
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
    spins_t, M, E = Ising.hybrid!(spins, β, 10^7)
    m = abs.(M[(length(M) ÷ 2):end]) / length(spins)
    mavg[k] = mean(m)
    mstd[k] = std(m)
end
scatter!(ax, βs, mavg, color=:red, markersize=5, label="MC, L=70")
errorbars!(ax, βs, mavg, mstd/2, color=:red, whiskerwidth=5)

axislegend(ax, position=:rb)
fig
```


## Wolff steps are mixed with Metropolis steps

```@example
import SquareIsingModel as Ising
using CairoMakie, Random

Random.seed!(1) # make reproducible
L = 50
β = Ising.βc
T = 20
Δ = 5
spins = Ising.random_configuration(L)
Ising.metropolis!(spins, β, 10^6) # equilibrate a bit to get clusters
spins_t, M, E = Ising.hybrid!(spins, β, T; save_interval = 1, local_steps = Δ)
fig = Figure(resolution=(600, 500))
for t ∈ 1:T
    ax = Axis(fig[cld(t, Δ), mod1(t, Δ)])
    hidedecorations!(ax)
    heatmap!(ax, spins_t[:,:,t], colormap=cgrad([:purple, :orange], [0.5]; categorical=true))
end
fig
```

In the figure above, each row is a sequence of consecutive Metropolis steps.
Inside a row the configurations are very similar differing in single sites.
When the row ends, a Wolff cluster move is taken.
It can be seen that the next row has suffered a larger change, since a cluster was flipped.


## Magnetic susceptibility

```@example
using Statistics, CairoMakie, Random
import SquareIsingModel as Ising

Random.seed!(1) # make reproducible
Ts = 2:0.01:3
βs = inv.(Ts)

fig = Figure(resolution=(600, 400))
ax = Axis(fig[1,1], xlabel=L"temperature $T$ ($=1/\beta$)", ylabel="susceptibility", yscale=log10)
@time for (L, color) in zip([4, 8, 16, 32], [:green, :orange, :blue, :red])
    χ = zeros(length(βs))
    for (k, β) in enumerate(βs)
        spins = Ising.random_configuration(L)
        spins_t, M, E = Ising.hybrid!(spins, β, 10^6)
        χ[k] = β/length(spins) * var(abs.(M))
    end
    scatter!(ax, Ts, χ, color=color, markersize=5, label=L"L=%$L")
end
vlines!(ax, [1 / Ising.βc], label=L"Onsager's $T_c$", color=:black, linestyle=:dash)
axislegend(ax, position=:rt)

fig
```


## Heat capacity

```@example
using Statistics, CairoMakie, Random
import SquareIsingModel as Ising

Random.seed!(1) # make reproducible
Ts = 1.5:0.01:3
βs = inv.(Ts)

fig = Figure(resolution=(600, 400))
ax = Axis(fig[1,1], xlabel=L"temperature $T$ ($=1/\beta$)", ylabel="heat capacity", yscale=log10)
@time for (L, color) in zip([4, 8, 16, 32], [:green, :orange, :blue, :red])
    C = zeros(length(βs))
    for (k, β) in enumerate(βs)
        spins = Ising.random_configuration(L)
        spins_t, M, E = Ising.hybrid!(spins, β, 10^7)
        C[k] = β^2/length(spins) * var(E)
    end
    scatter!(ax, Ts, C, color=color, markersize=5, label=L"L=%$L")
end
vlines!(ax, [1 / Ising.βc], label=L"Onsager's $T_c$", color=:black, linestyle=:dash)
axislegend(ax, position=:lt) 

fig
```


## Wolff vs. Metropolis spin flip rates

At low temperatures, Wolff flips more spins per unit time than Metropolis.
At high temperatures, Metropolis is more efficient.
The crossing point approaches the critical temperature for larger system sizes.

```@example
using Statistics, CairoMakie, Random
import SquareIsingModel as Ising

Random.seed!(1) # make reproducible
Ts = 1:0.5:4
βs = inv.(Ts)

fig = Figure(resolution=(600, 400))
ax = Axis(fig[1,1], xlabel=L"temperature $T$ ($=1/\beta$)", ylabel="spin flips / second", yscale=log10)
@time for (L, color) in zip([4, 8, 16, 32], [:blue, :red, :orange, :green])
    wolff_rates = zeros(length(βs))
    local_rates = zeros(length(βs))
    for (k, β) in enumerate(βs)
        spins = Ising.random_configuration(L)
        stats = Ising.HybridStats()
        spins_t, M, E = Ising.dynamic_hybrid!(spins, β, 10^6; hybrid_stats=stats)
        local_rates[k] = stats.local_flip / stats.local_time
        wolff_rates[k] = stats.wolff_flip / stats.wolff_time
    end
    scatterlines!(ax, Ts, local_rates, markersize=15, label=L"$L=%$L$ local", color=color, marker='x', markercolor=color, linestyle=:dash)
    scatterlines!(ax, Ts, wolff_rates, markersize=15, label=L"$L=%$L$ wolff", color=color, marker='o', markercolor=color, linestyle=:dot)
end
vlines!(ax, [1 / Ising.βc], label=L"Onsager's $T_c$", color=:black, linewidth=1)

Legend(fig[1,2], ax)

fig
```