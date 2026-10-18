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
    t, s = tanh(2β), sech(2β)^2
    jgK, _ = _onsager_elliptic(t, s)
    return -2t * (1 + 2s * jgK)
end

"""
    onsager_heat_capacity(β)

Specific heat capacity of the two-dimensional Ising model, derived by
Onsager in the thermodynamic limit.
"""
function onsager_heat_capacity(β::Real)
    t, s = tanh(2β), sech(2β)^2
    isinf(β) && return zero(t) # zero temperature, where β * s would be Inf * 0
    jgK, h = _onsager_elliptic(t, s)
    return 4(β * s)^2 * (4t^2 * h - 2jgK)
end

#=
Onsager's internal energy and heat capacity are usually written as

    u = -coth(2β) * (1 + 2/π * j * K)
    c = β^2 * coth(2β)^2 * (2/π) * (((j - 1/2)^2 + 7/4) * K - 2E - (1 - j) * π/2)

where j = 2tanh(2β)^2 - 1, and K, E are the complete elliptic integrals of parameter
m = k^2, k = 2tanh(2β) / cosh(2β). Evaluated as written, these give 0 * Inf at β = 0, and lose
precision to cancellations at high temperature (and c also at low temperature), where m → 0
and 2/π * K → 1, 2/π * E → 1. With t = tanh(2β), s = sech(2β)^2, so that j = t^2 - s and m = 4t^2 s,
they are rewritten above as

    u = -2t * (1 + 2s * j * gK)
    c = 4β^2 s^2 * (4t^2 * h - 2j * gK)

in terms of gK = (2/π * K - 1) / m and h = (4/π * (K - E) - m) / m^2, which are smooth at m = 0
and are computed from their Taylor series when m is small. Near βc, m → 1 and gK diverges,
but j * gK → 0.
=#
const _ONSAGER_K_COEFS = ntuple(n -> (binomial(2n, n) / 4^n)^2, 21) # 2/π * K = 1 + Σ_n K_n m^n
const _ONSAGER_GK_COEFS = _ONSAGER_K_COEFS[1:20]
const _ONSAGER_H_COEFS = ntuple(n -> 2 * _ONSAGER_K_COEFS[n + 1] * (2n + 2) / (2n + 1), 20)
function _onsager_elliptic(t::Real, s::Real)
    j = t^2 - s
    m = min(4t^2 * s, 1) # may otherwise round above 1 near βc
    if m < 0.1
        gK = evalpoly(m, _ONSAGER_GK_COEFS)
        h = evalpoly(m, _ONSAGER_H_COEFS)
    else
        K = SpecialFunctions.ellipk(m)
        if isinf(K) # near βc, m rounds to 1, where K ≈ log(4/|j|)
            K = log(4 / abs(j))
        end
        E = SpecialFunctions.ellipe(m)
        gK = (2/π * K - 1) / m
        h = (4/π * (K - E) - m) / m^2
    end
    jgK = iszero(j) ? zero(j) : j * gK
    return jgK, h
end

"""
    kramers_wannier(β)

Returns the Kramers-Wannier dual inverse-temperature of β.
"""
kramers_wannier(β::Real) = -log(tanh(β)) / 2
