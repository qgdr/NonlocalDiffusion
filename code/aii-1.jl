using Plots
using OffsetArrays
using LinearAlgebra



# f(x, y, α) = x*y *(1+ x^(2-α) + y^(2-α) ) + (1+x+y)^(3-α) + (1+x)*(1+y) * (1- (1+x)^(2-α) - (1+y)^(2-α) )


# f(8//10, 14//10, 1.01)

# plot(y->f(10//10, y,  1.5), 1, 1000)
# plot(x->f(x, 101//100,  1.1), 0, 10)

# df_dx(x, y, α) = 1 + 2y + y^(3-α) - (1+y)^(3-α)  +  (3-α)* ( x^(2-α) * y + (1+x+y)^(2-α) -(1+x)^(2-α)*(1+y) )
# plot(x->df_dx(x, 12//10,  1.01), 0, 10)

# plot(y->df_dx(1, y,  1.01), 0, 2)

########################

function delta(x)
    # a = Ω[1]
    # b = Ω[2]
    return min(abs(x), abs(1-x))
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

Ω = (0, 1)


α = 1+ 1//10 #19999 // 10000
α = BigFloat(α)

@show 1 < α < 2

r = 1 # 4/α
N = 100
x = Gridmesh(Ω, N, r)        # [x0, x1, ..., x2N]
h = x |> parent |> diff      # [h1, ...h2N], h_j = x_j - x_{j-1}

s(i) = (abs(x[i] - x[1])^(3-α) - abs(x[i] - x[0])^(3-α) ) / h[1]
p(i) = 2 .* abs(1/2 .- x[i])^(3-α)
q(i) = x[i]^(3-α) 

Si = zeros(BigFloat, 2N-1)
for i=1:2N-1
    Si[i] = (s(i+1) / h[i+1] / (h[i] + h[i+1]) - s(i) / (h[i] * h[i+1]) + s(i-1) / h[i] / (h[i] + h[i+1]))
end

Pi = zeros(BigFloat, 2N-1)
for i=1:2N-1
    Pi[i] = (p(i+1) / h[i+1] / (h[i] + h[i+1]) - p(i) / (h[i] * h[i+1]) + p(i-1) / h[i] / (h[i] + h[i+1]))
end

Qi = zeros(BigFloat, 2N-1)
for i=1:2N-1
    Qi[i] = (q(i+1) / h[i+1] / (h[i] + h[i+1]) - q(i) / (h[i] * h[i+1]) + q(i-1) / h[i] / (h[i] + h[i+1]))
end

# plot(Si)
plot(Si .* x[1:2N-1].^α ./ ( (3-α)*(2-α)*(α-1 ) ))
# plot(Pi)
plot!(Pi .* (1//2 .- delta.(x[1:2N-1]) .+ 1//2N).^(α-1) ./ ( (3-α)*(2-α) ) )
# plot!(Qi)
plot!(Qi .* x[1:2N-1].^(α-1) ./ ( (3-α)*(2-α) )  )
ylims!(0, 1)

# λ = 1//100
# M = λ * Si + Pi - Qi
# @show all(M.> 0)
# plot(x[1:2N-1], M)

