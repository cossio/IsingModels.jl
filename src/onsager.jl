"""
    βc

Critical temperature of the 2-dimensional infinite lattice Ising model determined
by Onsager.
"""
const βc = log1p(√2) / 2

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
    j = 2tanh(2β)^2 - 1
    m = 1 - j^2 # = k^2, k = 2tanh(2β) / cosh(2β), but unlike k^2 can't round above 1
    K = SpecialFunctions.ellipk(m)
    if isinf(K)
        # near βc, k^2 rounds to 1 and K ≈ log(4/|j|) diverges, but j * K → 0
        jK = iszero(j) ? zero(j) : j * log(4 / abs(j))
    else
        jK = j * K
    end
    return -coth(2β) * (1 + 2/π * jK)
end

"""
    onsager_heat_capacity(β)

Specific heat capacity of the two-dimensional Ising model, derived by
Onsager in the thermodynamic limit.
"""
function onsager_heat_capacity(β::Real)
    j = 2tanh(2β)^2 - 1
    m = 1 - j^2 # = k^2, k = 2tanh(2β) / cosh(2β), but unlike k^2 can't round above 1
    K = SpecialFunctions.ellipk(m)
    E = SpecialFunctions.ellipe(m)
    return β^2 * coth(2β)^2 * (2/π) * (((j - 1//2)^2 + 7//4) * K - 2E - (1 - j) * π / 2)
end

"""
    kramers_wannier(β)

Returns the Kramers-Wannier dual inverse-temperature of β.
"""
kramers_wannier(β::Real) = -log(tanh(β)) / 2
