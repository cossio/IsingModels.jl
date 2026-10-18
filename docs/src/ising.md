# The 2-dimensional Ising model

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

In the thermodynamic limit (large `L` with `K = L`), this model suffers a phase transition at the critical inverse temperature $\beta \approx 0.44$ (called [`βc`](@ref Ising.βc) in the package).

In this package, the system is simulated using the Metropolis algorithm or the Wolff cluster algorithm, both explained here:

```
Newman, Mark EJ, and G. T. Barkema. "Monte Carlo Methods in Statistical Physics (1999)." New York: Oxford 475 (1999).
```

Onsager derived exact expressions for the free energy, the heat capacity, and the internal energy in the thermodynamic limit.

```
Onsager, Lars. "Crystal statistics. I. A two-dimensional model with an order-disorder transition." Physical Review 65.3-4 (1944): 117.
```

See [`onsager_internal_energy`](@ref Ising.onsager_internal_energy) and [`onsager_heat_capacity`](@ref Ising.onsager_heat_capacity), implemented in this package.

## Spin configurations

Spin configurations are arrays of `Bool`, where `true` is spin ``+1`` and `false` is spin ``-1`` (see [`spin`](@ref Ising.spin) and [`binary`](@ref Ising.binary)).
Any `AbstractArray{Bool}` works: a `BitArray` (from `bitrand`, `trues` or `falses`), an `Array{Bool}` (from `rand(Bool, L, K)`), or a view of either.

A `BitArray` stores one spin per bit, using 8 times less memory than an `Array{Bool}`, which stores one spin per byte.
Reading or writing a single spin takes a bit more work in a `BitArray`, so an `Array{Bool}` is usually faster for lattices that fit in the CPU cache.
For very large lattices, the smaller memory footprint of a `BitArray` can make it faster instead.
The configurations saved during a simulation are always stored in a `BitArray`.
