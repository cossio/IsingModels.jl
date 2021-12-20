# Examples using the Metropolis+F(M) sampling method

#=
The `metropolis!(...; f = f)` method samples an energy of the form:

```math
E = -\sum_{\langle i,j\rangle} s_i s_j + f(M)
```

where ``f(M)`` is a function of the total magnetization,

```math
M = \sum_i s_i
```

Note that ``M`` is not normalized by the number of spins.
Also note that the temperature multiplies both the original Ising energy, and ``f(M)``:

```math
- \beta E = \beta\sum_{\langle i,j\rangle} s_i s_j - \beta f(M)
```

Therefore the system prefers configurations with *lower* values of `f(M)`.

We will consider the above energy with the modified term ``f(M) = w|M|/\beta``.
Here ``|M|`` is the absolute value of the magnetization.
``w`` is a factor that weights this term in the energy, and we divide by ``\beta``,
so that the overall system looks like this:

```math
- \beta E = \beta\sum_{\langle i,j\rangle} s_i s_j - w|M|
```

First load required packages.
=#

import IsingModels as Ising
using Statistics, Random
using LogExpFunctions, CairoMakie, IrrationalConstants
nothing #hide

# Try to make output reproducible

Random.seed!(1)
nothing #hide

# Define the parameter ranges we will consider.

ws = 0:0.001:0.050 # weight of extra term
Ls = [32, 64]
βs = [0.4, 0.5, 0.6] # inverse temperatures
cs = [:red, :purple, :blue] # colors
nothing #hide

# Simulate and collect data.

magnetization_avg = Dict{typeof((β=first(βs), w=first(ws), L=first(Ls))), Float64}()
magnetization_std = Dict{typeof((β=first(βs), w=first(ws), L=first(Ls))), Float64}()

for β in βs, w in ws, L in Ls
    f(M::Real) = w * abs(M) / β
    σ = bitrand(L, L)
    σ_t, M, E = Ising.metropolis!(σ, β; steps=10^7, f=f)
    magnetization_avg[(β=β, w=w, L=L)] = mean(M)
    magnetization_std[(β=β, w=w, L=L)] = std(M)
end

nothing #hide

# Now let's plot the results.

fig = Figure(resolution=(1000, 400))
for (iL, L) in enumerate(Ls)
    ax = Axis(fig[1,iL], xlabel="w", ylabel="m", title="L=$L")
    for (iβ, (β, color)) in enumerate(zip(βs, cs))
        mavg = [magnetization_avg[(β=β, w=w, L=L)] / L^2 for w in ws]
        mstd = [magnetization_std[(β=β, w=w, L=L)] / L^2 for w in ws]
        lines!(ax, ws, mavg, label="β=$β", color=color)
        errorbars!(ax, ws, mavg, mstd/2, whiskerwidth=5, color=color)
    end
    axislegend(ax, position=:rt)
end
fig

#=
We now try the function ``f(M) = \log\cosh(wM) / \beta``, so that:

```math
- \beta E = \beta\sum_{\langle i,j\rangle} s_i s_j - \log\cosh(w M)
```
=#

# First collect some data.

magnetization_avg = Dict{typeof((β=first(βs), w=first(ws), L=first(Ls))), Float64}()
magnetization_std = Dict{typeof((β=first(βs), w=first(ws), L=first(Ls))), Float64}()

for β in βs, w in ws, L in Ls
    f(M) = logcosh(w * M) / β
    σ = bitrand(L, L)
    σ_t, M, E = Ising.metropolis!(σ, β; steps=10^7, f=f)
    magnetization_avg[(β=β, w=w, L=L)] = mean(M)
    magnetization_std[(β=β, w=w, L=L)] = std(M)
end

nothing #hide

# Now plot the result.

fig = Figure(resolution=(1000, 400))
for (iL, L) in enumerate(Ls)
    ax = Axis(fig[1,iL], xlabel="w", ylabel="m", title="L=$L")
    for (iβ, (β, color)) in enumerate(zip(βs, cs))
        mavg = [magnetization_avg[(β=β, w=w, L=L)] / L^2 for w in ws]
        mstd = [magnetization_std[(β=β, w=w, L=L)] / L^2 for w in ws]
        lines!(ax, ws, mavg, label="β=$β", color=color)
        errorbars!(ax, ws, mavg, mstd/2, whiskerwidth=5, color=color)
    end
    axislegend(ax, position=:rt)
end
fig
