"""
    βc

Critical temperature of the 2-dimensional infinite lattice Ising model determined
by Onsager.
"""
const βc = log(1 + √big(2)) / 2

"""
    onsager_magnetization(β)

Analytical magnetization of the two-dimensional Ising model
found by Onsager in the thermodynamic limit.
"""
function onsager_magnetization(β::Real)
    @assert β ≥ 0
    return max(1 - csch(2β)^4, 0)^(1/oftype(β, 8))
end

"""
    onsager_internal_energy(β)

Internal energy per site of the two-dimensional Ising model, derived by
Onsager in the thermodynamic limit.
"""
function onsager_internal_energy(β::Real)
    k = 2sinh(2β) / cosh(2β)^2
    return -2coth(2β) * (1 + 2/π * (2tanh(2β)^2 - 1) * ellipk(k^2))
end

"""
    onsager_heat_capacity(β)

Specific heat capacity of the two-dimensional Ising model, derived by
Onsager in the thermodynamic limit.
"""
function onsager_heat_capacity(β::Real)
    k = 2sinh(2β) / cosh(2β)^2
    k_ = 2tanh(2β)^2 - 1
    β * coth(2β)^2 * 2/π * (2ellipk(k^2) - 2ellipe(k^2) - (1 - k_) * (π/2 + k_ * ellipk(k^2)))
end

"""
    kramers_wannier(β)

Returns the Kramers-Wannier dual inverse-temperature of β.
"""
kramers_wannier(β::Real) = -log(tanh(β)) / 2
