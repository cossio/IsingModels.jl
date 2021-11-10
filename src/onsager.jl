"""
    βc

Critical temperature of the 2-dimensional infinite lattice Ising model determined
by Onsager.
"""
const βc = log(1 + √big(2)) / 2

"""
    onsager_magnetization(β)

Analytical expression for the magnetization found by Onsager in the thermodynamic limit.
"""
function onsager_magnetization(β::Real)
    @assert β ≥ 0
    return max(1 - csch(β)^4, 0)^(1/oftype(β, 8))
end
