# Reference

This package doesn't export any symbols. 

## Grid lattice

```@docs
SquareIsingModel.energy
SquareIsingModel.adjacency_matrix
SquareIsingModel.neighbors
SquareIsingModel.neighbor_sum
SquareIsingModel.random_configuration
```

## Onsager analytical solution

```@docs
SquareIsingModel.βc
SquareIsingModel.onsager_magnetization
SquareIsingModel.onsager_internal_energy
SquareIsingModel.onsager_heat_capacity
SquareIsingModel.kramers_wannier
```

## Monte Carlo simulations

```@docs
SquareIsingModel.metropolis!
SquareIsingModel.wolff!
SquareIsingModel.hybrid!
```