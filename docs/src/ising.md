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

In the thermodynamic limit (large `L` with `K = L`), this model suffers a phase transition at the critical inverse temperature $\beta \approx 0.44$ (called [`βc`](@ref) in the package).
