using LinearAlgebra
using OffsetArrays

using Plots



r = 300 // 100
Ω = (0, 1)
α = 15 // 10
α = BigFloat(α)

@show 1 < α < 2

function Dh2T(N, i, j, θ)
    # global α

    y(j, θ) = θ * (j - 1)^r + (1 - θ) * j^r
    h(j) = j^r - (j - 1)^r

    # yjmθ = θ * (j-2)^r + (1 - θ) * (j-1)^r
    # yjpθ = θ * j^r + (1 - θ) * (j+1)^r


    K(i, j) = h(j)^3 * abs(i^r - y(j, θ))^(1 - α) * y(j, θ)^(α / 2 - 2)

    # Kjm = h[j-1]^3 * abs(x[i-1] - yjmθ)^(1 - α) * yjmθ^(α / 2 - 2)
    # Kjp = h[j+1]^3 * abs(x[i+1] - yjpθ)^(1 - α) * yjpθ^(α / 2 - 2)

    return 2 / (h(i) + h(i + 1)) * (1 / h(i + 1) * K(i + 1, j + 1) - (1 / h(i + 1) + 1 / h(i)) * K(i, j) + 1 / h(i) * K(i - 1, j - 1))
end

Ns = 2 .^ (6:18)

result = []

for N ∈ Ns


    i = N ÷ 2
    # j = i ÷ 2
    # j = i
    j = i - 2
    # j = 3i ÷ 4

    push!(result, Dh2T(N, i, j, 0.4))

end





# result = result .* Ns.^2
# diff( abs.(result) .|> log)./diff(log.(Ns))


plot(log.(Ns), log.(abs.(result)))

# plot(result)
