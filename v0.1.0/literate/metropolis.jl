# Examples using the Metropolis sampling method

## Magnetization as a function of temperature

import IsingModels as Ising
using Statistics, CairoMakie, Random

Random.seed!(1) # make reproducible
nothing #hide

# Define the temperatures we will consider.

βs = 0:0.025:1
nothing #hide

#=
Simulate a range of temperatures.
Compare to analytical magnetization.
=#

fig = Figure(resolution=(600,400))
ax = Axis(fig[1,1], xlabel="β", ylabel="m")
lines!(ax, 0:0.01:1, Ising.onsager_magnetization, color=:black, label="analytical")

mavg = zeros(length(βs))
mstd = zeros(length(βs))
σ = bitrand(20, 20)
for (k, β) in enumerate(βs)
    σ_t, M, E = Ising.metropolis!(σ, β; steps=10^7)
    m = abs.(M[(length(M) ÷ 2):end]) / length(σ)
    mavg[k] = mean(m)
    mstd[k] = std(m)
end
scatter!(ax, βs, mavg, color=:blue, markersize=5, label="MC, L=20")
errorbars!(ax, βs, mavg, mstd/2, color=:blue, whiskerwidth=5)

mavg = zeros(length(βs))
mstd = zeros(length(βs))
σ = bitrand(64, 64)
for (k, β) in enumerate(βs)
    σ_t, M, E = Ising.metropolis!(σ, β; steps=10^7)
    m = abs.(M[(length(M) ÷ 2):end]) / length(σ)
    mavg[k] = mean(m)
    mstd[k] = std(m)
end
scatter!(ax, βs, mavg, color=:red, markersize=5, label="MC, L=50")
errorbars!(ax, βs, mavg, mstd/2, color=:red, whiskerwidth=5)

axislegend(ax, position=:rb)
fig
