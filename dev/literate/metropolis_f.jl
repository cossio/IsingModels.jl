# Examples using the Metropolis+F(M) sampling method

#=
The `metropolis_f!` method samples an energy of the form:

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

Let's look at some examples.

First load required packages.
=#

import SquareIsingModel as Ising
using Statistics, Random
using LogExpFunctions, CairoMakie, IrrationalConstants
nothing #hide

# Try to make output reproducible

Random.seed!(1)
nothing #hide

# Define the temperatures we will consider.

βs = 0:0.025:1
nothing #hide

#=
We will consider the above energy with the modified term ``f(M) = w|M|/\beta``.
Here ``|M|`` is the absolute value of the magnetization.
``w`` is a factor that weights this term in the energy, and we divide by ``\beta``,
so that the overall system looks like this:

```math
- \beta E = \beta\sum_{\langle i,j\rangle} s_i s_j - w|M|
```
=#

#=
We now simulate a range of values of `w` and some system sizes.
=#

ws = [0 1e-2; 1e-1 1; 10 100] # values of w
Ls = [32, 64] # system size
nothing #hide

# Simulate!

colors = [:blue, :red]
fig = Figure(resolution=(1000, 150 * length(ws)))
for iw in CartesianIndices(ws),
    w = ws[iw]
    ax = Axis(fig[Tuple(iw)...], xlabel="β", ylabel="m", title="w=$w")
    lines!(ax, 0:0.01:1, Ising.onsager_magnetization, color=:black, label="Onsager's M")

    for (iL, L) in enumerate(Ls)
        mavg = zeros(length(βs))
        mstd = zeros(length(βs))
        spins = Ising.random_configuration(L)
        for (k, β) in enumerate(βs)
            f(M) = w * abs(M) / β
            spins_t, M, E = Ising.metropolis_f!(spins, β; steps=10^7, f=f)
            m = abs.(M[(length(M) ÷ 2):end]) / length(spins)
            mavg[k] = mean(m)
            mstd[k] = std(m)
        end
        scatter!(ax, βs, mavg, color=colors[iL], markersize=5, label="L=$L")
        errorbars!(ax, βs, mavg, mstd/2, color=colors[iL], whiskerwidth=5)
    end

    if Tuple(iw) == (1,1)
        axislegend(ax, position=:lt)
    end
end
fig

#=
We now try the function ``f(M) = \log\cosh(wM) / \beta``, so that:

```math
- \beta E = \beta\sum_{\langle i,j\rangle} s_i s_j - \log\cosh(w M)
```
=#

#=
First let's define a function to compute `log(cosh(x))` in a numerically stable way.
=#

function logcosh(x::Real)
    abs_x = abs(x)
    return abs_x + log1pexp(-2 * abs_x) - logtwo
end

# Now we are ready to run the simulation.

colors = [:blue, :red]
fig = Figure(resolution=(1000, 150 * length(ws)))
for iw in CartesianIndices(ws)
    w = ws[iw]
    ax = Axis(fig[Tuple(iw)...], xlabel="β", ylabel="m", title="w=$w")
    lines!(ax, 0:0.01:1, Ising.onsager_magnetization, color=:black, label="Onsager's M")

    for (iL, L) in enumerate(Ls)
        mavg = zeros(length(βs))
        mstd = zeros(length(βs))
        spins = Ising.random_configuration(L)
        for (k, β) in enumerate(βs)
            f(M) = logcosh(w * M) / β
            spins_t, M, E = Ising.metropolis_f!(spins, β; steps=10^7, f=f)
            m = abs.(M[(length(M) ÷ 2):end]) / length(spins)
            mavg[k] = mean(m)
            mstd[k] = std(m)
        end
        scatter!(ax, βs, mavg, color=colors[iL], markersize=5, label="L=$L")
        errorbars!(ax, βs, mavg, mstd/2, color=colors[iL], whiskerwidth=5)
    end

    if Tuple(iw) == (1,1)
        axislegend(ax, position=:lt)
    end
end
fig
