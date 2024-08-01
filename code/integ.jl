# using SymPy

# x, y, α = symbols("x, y, alpha")
# @show integrate(v, (y, 0, 1))


using Cubature
using SpecialFunctions

α = BigFloat(11//10)
x = BigFloat(0.2)


# function v(y)
#     return ( y * (1-y) )^(α/2) * abs(x - y)^ (1 - α)
# end


# hquadrature(x -> x^3, 0, 1)


function u_exact(x, α)
    α = BigFloat(α)
    # CC = 2^(-α) * gamma(1 / 2) / (gamma(1 + α / 2) * gamma((1 + α) / 2))
    CC=1
    ut = CC * ( x * (1 - x) ) ^ (α / 2)
    return ut
end

function v(y)
    return u_exact(y, α) * abs(x - y)^(1-α)
end

@show (val, err) = hquadrature(v, 0, 1, reltol=1e-12);

using OffsetArrays
b = OffsetArray([3, 4 ,5], -1:1)
for i in eachindex(b)
    println(i)
end
