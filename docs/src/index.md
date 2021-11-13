# SquareIsingModel.jl Documentation

## The 2-dimensional Ising model

The 2-dimensional Ising model is defined by the energy function:

```math
E(\mathbf{\sigma}) = - \sum_{\langle i j \rangle} \sigma_i \sigma_j
```

where $\langle i j \rangle$ refers to connected pairs of sites in the square grid lattice, and $\sigma_i = \pm 1$ are spins.
At inverse temperature $\beta$, this defines a Boltzmann probability distribution:

```math
P(\mathbf{\sigma}) = \frac{1}{Z} \mathrm{e}^{-\beta E (\mathbf{\sigma})}
```

where

```math
Z = \sum_{\mathbf{\sigma}} \mathrm{e}^{-\beta E(\mathbf{\sigma})}
```

is the partition function.

In the two-dimensional grid lattice, we assume we have a $L\times K$ plane grid, where each spin is connected to its four neighbors.
We assume periodic boundary conditions, so spin `(1,1)` is connected to spin `(L,K)`.

In the thermodynamic limit (large `L` with `K = L`), this model suffers a phase transition at the critical inverse temperature $\beta \approx 0.44$.

```@example
using SquareIsingModel: metropolis!, onsager_magnetization, βc
using Statistics, CairoMakie, Random

Random.seed!(1) # make reproducible

βs = 0:0.05:1
mavg = zeros(length(βs))
mstd = zeros(length(βs))
spins = rand((Int8(1), Int8(-1)), 50, 50);
for (k,β) in enumerate(βs)
    println("simulating β = $β ...")
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

## Reference

This package doesn't export any symbols. 

### Grid lattice

```@docs
SquareIsingModel.energy
SquareIsingModel.adjacency_matrix
SquareIsingModel.neighbors
SquareIsingModel.neighbor_sum
SquareIsingModel.random_configuration
```

### Onsager analytical solution

```@docs
SquareIsingModel.βc
SquareIsingModel.onsager_magnetization
```

### Monte Carlo simulations

```@docs
SquareIsingModel.metropolis!
SquareIsingModel.wolff!
```
