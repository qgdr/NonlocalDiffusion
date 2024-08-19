using SpecialFunctions
using LinearAlgebra
# using GR
using XLSX
using OffsetArrays
using Cubature
using Plots


function function_f(x, α)
    f = 1
    return f
end

function u_exact(x, α)
    α = BigFloat(α)
    CC = 2^(-α) * gamma(1 / 2) / (gamma(1 + α / 2) * gamma((1 + α) / 2))
    ut = CC * (x * (1 - x))^(α / 2)
    return ut
end


# function Gridmesh(Ω, N, r)
#     grid = OffsetArray(zeros(BigFloat, 2N + 1), 0:2N)
#     a = Ω[1]
#     b = Ω[2]

#     for j = 0:N
#         grid[j] = (BigFloat((b - a)) / 2) * (j // N)^r .+ a
#     end
#     for j = N+1:2N
#         grid[j] = -(BigFloat((b - a)) / 2) * (2 - j // N)^r .+ b
#     end
#     # end

#     return grid
# end





r = 15 // 10
Ω = (0, 1)
a = Ω[1]
b = Ω[2]

α = 15 // 10
α = BigFloat(α)

@show 1 < α < 2

R1 = []
RN = []


Ns = 2 .^(10:15)
for N ∈ Ns

    function x(j)
        if j ∈ 0:N
            return (BigFloat((b - a)) / 2) * (j // N)^r .+ a
        else
            # for j = N+1:2N
            return -(BigFloat((b - a)) / 2) * (2 - j // N)^r .+ b
        end
    end
    h(j) = x(j) - x(j - 1)


    CR = 1 / (2 * cos(π * α / 2) * gamma(2 - α))

    ## D_{ij} = |x_i-x_j|^{3-α}



    D(i, j) = abs(x(i) - x(j))^(3 - α)


    ## ̃a_{ij} = Cₐ (|x_{j+1} - x_i|^{3-α} / h_{j+1} - (h_j+h_{j+1})/(h_j*h_{j+1}) * |x_j - x_i|^{3-α} +  |x_{j-1} - x_i|^{3-α}/h_j)


    C_a = 1 / (2 - α) / (3 - α)

    A_1(i, j) = C_a * (D(i, j + 1) / h(j + 1) - (h(j) + h(j + 1)) / (h(j) * h(j + 1)) * D(i, j) + D(i, j - 1) / h(j))


    ## a_{ij} = 2 Cᵣ ( ̃a_{i+1,j}/ h_{i+1} - (h_i+h_{i+1})/(h_i*h_{i+1}) ̃a_{ij} +  a_{i-1,j}/ h_{i})




    A(i, j) = 2 * CR * (A_1(i + 1, j) / h(i + 1) / (h(i) + h(i + 1)) - A_1(i, j) / (h(i) * h(i + 1)) + A_1(i - 1, j) / h(i) / (h(i) + h(i + 1)))




    # truncation error

    U(j) = u_exact.(x(j), α)

    # i = 4N ÷ 5
    # i = 2
    i= N-3

    FN = 0
    for j = 1:2N-1
        FN += A(i, j) * U(j)
    end

    R = FN - 1
    push!(RN, R)


end

plot(log.(Ns), abs.(RN) .|> log)
#  diff( abs.(RN) .|> log)./diff(log.(Ns))

