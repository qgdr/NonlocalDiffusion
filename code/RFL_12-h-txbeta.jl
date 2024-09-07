using SpecialFunctions
using LinearAlgebra
# using GR
using XLSX
using OffsetArrays
using Cubature
using Plots

β = -18//10
β = BigFloat(β)
function f(x, α)
    return x^β * exp(x)
end

# function u_exact(x, α)
#     α = BigFloat(α)
#     CC = 2^(-α) * gamma(1 / 2) / (gamma(1 + α / 2) * gamma((1 + α) / 2))
#     ut = CC * (x * (1 - x))^(α / 2)
#     return ut
# end


function Gridmesh(Ω, N, r)
    grid = OffsetArray(zeros(BigFloat, 2N + 1), 0:2N)
    a = Ω[1]
    b = Ω[2]

    for j = 0:N
        grid[j] = (BigFloat((b - a)) / 2) * (j // N)^r .+ a
    end
    for j = N+1:2N
        grid[j] = -(BigFloat((b - a)) / 2) * (2 - j // N)^r .+ b
    end
    # end

    return grid
end





r = 30 // 10
Ω = (0, 1)


α = 19 // 10
α = BigFloat(α)

@show 1 < α < 2

UN = []



for N ∈ [100, 200, 400]

# function num_err(Ω, r, α, N)

x = Gridmesh(Ω, N, r)        # [x0, x1, ..., x2N]
h = x |> parent |> diff      # [h1, ...h2N], h_j = x_j - x_{j-1}





CR = 1 / (2 * cos(π * α / 2) * gamma(2 - α))

## D_{ij} = |x_i-x_j|^{3-α}

D = zeros(BigFloat, 2N + 1, 2N + 1)
D = OffsetArray(D, 0:2N, 0:2N)

for i = 0:2N, j = 0:2N
    D[i, j] = abs(x[i] - x[j])^(3 - α)
end

## ̃a_{ij} = Cₐ (|x_{j+1} - x_i|^{3-α} / h_{j+1} - (h_j+h_{j+1})/(h_j*h_{j+1}) * |x_j - x_i|^{3-α} +  |x_{j-1} - x_i|^{3-α}/h_j)

A_1 = OffsetArray(zeros(BigFloat, 2N + 1, 2N - 1), 0:2N, 1:2N-1)

C_a = 1 / (2 - α) / (3 - α)
for j = 1:2N-1
    A_1[:, j] = C_a * (D[:, j+1] / h[j+1] - (h[j] + h[j+1]) / (h[j] * h[j+1]) * D[:, j] + D[:, j-1] / h[j])
end

## a_{ij} = 2 Cᵣ ( ̃a_{i+1,j}/ h_{i+1} - (h_i+h_{i+1})/(h_i*h_{i+1}) ̃a_{ij} +  a_{i-1,j}/ h_{i})

A = zeros(BigFloat, 2N - 1, 2N - 1)

for i = 1:2N-1
    A[i, :] = 2 * CR * (A_1[i+1, :] / h[i+1] / (h[i] + h[i+1]) - A_1[i, :] / (h[i] * h[i+1]) + A_1[i-1, :] / h[i] / (h[i] + h[i+1]))
end


U_s = parent(A) \ f.(x[1:2N-1], α)

push!(UN, U_s)


end

# UN
# dU1 = zeros(BigFloat, 50)
# dU2 = zeros(BigFloat, 50)

dU1 = maximum(abs.(UN[2][2:2:end] - UN[1]))
dU2 = maximum(abs.(UN[3][4:4:end] - UN[2][2:2:end]))

log2.(dU1 ./ dU2)


# (abs(UN[2][100]-UN[1][50]) / abs(UN[3][200]-UN[2][100])) |> log2