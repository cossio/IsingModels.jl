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
    j = 2tanh(2β)^2 - 1
    K = ellipk(k^2)
    return -coth(2β) * (1 + 2/π * j * K)
end

"""
    onsager_heat_capacity(β)

Specific heat capacity of the two-dimensional Ising model, derived by
Onsager in the thermodynamic limit.
"""
function onsager_heat_capacity(β::Real)
    k = 2sinh(2β) / cosh(2β)^2
    K = ellipk(k^2)
    E = ellipe(k^2)
    j = 2tanh(2β)^2 - 1
    β^2 * coth(2β)^2 * (2/π) * (2K - 2E - (1 - j) * (π/2 + j * K))
end

"""
    kramers_wannier(β)

Returns the Kramers-Wannier dual inverse-temperature of β.
"""
kramers_wannier(β::Real) = -log(tanh(β)) / 2
