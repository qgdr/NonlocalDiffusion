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





r = 11//10
Ω = (0, 1)


α = 11 // 10
α = BigFloat(α)

@show 1 < α < 2

R1 = []
RN = []



for N ∈ 100:100:1000


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




# truncation error

U = u_exact.(x, α)
F = parent(A) * U[1:2N-1]

# plot(x[1:2N-1], F)

# U_s = parent(A) \ ones(2N - 1)
# plot(x[1:2N-1], U[1:2N-1]-U_s)
# plot!(x[1:2N-1], U_s)

Si = sum(parent(A), dims=2)









trunc_err = abs.(F .- 1) |> parent |> maximum
# # u_err = abs.(U[1:2N-1] - U_s) |> maximum
# R = parent(abs.(F .- 1)) .* (x[1:2N-1] .^ α)
# # println(Ri[l])
# Ri = maximum(R[1:N])
# S = maximum( abs.(Si[1:N] .* x[1:N] .* α) )



# plot(trunc_err)
# plot(Si ./ (x[1:2N-1] .^ (-α) .+ (1 .- x[1:2N-1]) .^ (-α)) )
# plot(F .- 1)

## tunc_err
R = F .- 1
tar_odr = R ./ Si

# tar_odr[1]
push!(R1, R[1])
push!(RN, R[N-3])


end

plot(log.(1:10), abs.(RN) .|> log)
#  diff( abs.(RN) .|> log)./diff(log.(1:10))

