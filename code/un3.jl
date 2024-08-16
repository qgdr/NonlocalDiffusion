using LinearAlgebra
using OffsetArrays

using Plots


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


r = 11 // 10
Ω = (0, 1)
α = 13 // 10
α = BigFloat(α)

@show 1 < α < 2

function Dh2T(i, j, θ, x, h)
    global α

    yjθ = θ * x[j-1] + (1 - θ) * x[j]

    yjmθ = θ * x[j-2] + (1 - θ) * x[j-1]
    yjpθ = θ * x[j] + (1 - θ) * x[j+1]

    Kj = h[j]^3 * abs(x[i] - yjθ)^(1 - α) * yjθ^(α / 2 - 2)

    Kjm = h[j-1]^3 * abs(x[i-1] - yjmθ)^(1 - α) * yjmθ^(α / 2 - 2)
    Kjp = h[j+1]^3 * abs(x[i+1] - yjpθ)^(1 - α) * yjpθ^(α / 2 - 2)

    return 2 / (h[i] + h[i+1]) * (1 / h[i+1] * Kjp - (1 / h[i+1] + 1 / h[i]) * Kj + 1 / h[i] * Kjm)
end

Ns = 2 .^ (6:17)

result = []

for N ∈ Ns

    x = Gridmesh(Ω, N, r)        # [x0, x1, ..., x2N]
    h = x |> parent |> diff      # [h1, ...h2N], h_j = x_j - x_{j-1}

    i = N - 2
    j = N - 2

    push!(result, Dh2T(i, j, 0.5, x, h))

end






# diff( abs.(result) .|> log)./diff(log.(Ns))


plot(log.(Ns), log.( abs.(result)))

# plot(result)
