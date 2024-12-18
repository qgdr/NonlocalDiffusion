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





Ω = (0, 1)


@show α = 18 // 10
α = BigFloat(α)
@show 1 < α < 2

@show r = 4/α


function num_err(Ω, r, α, N)

    x = Gridmesh(Ω, N, r)        # [x0, x1, ..., x2N]
    h = x |> parent |> diff      # [h1, ...h2N], h_j = x_j - x_{j-1}



    ## \int_Ω u(y) / |x-y|^{α-1} dy

    # int_uoa, int_err = OffsetArray(zeros(BigFloat, 2N+1), 0:2N), OffsetArray(zeros(BigFloat, 2N+1), 0:2N)

    # Threads.@threads for i ∈ eachindex(x)
    #     v(y) = u_exact(y, α) * abs(x[i] - y)^(1-α)
    #     (val, err) = hquadrature(v, 0, 1, reltol=1e-12)
    #     int_uoa[i], int_err[i] = val, err    
    # end

    # A

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
    Threads.@threads for j = 1:2N-1
        A_1[:, j] = C_a * (D[:, j+1] / h[j+1] - (h[j] + h[j+1]) / (h[j] * h[j+1]) * D[:, j] + D[:, j-1] / h[j])
    end

    ## a_{ij} = 2 Cᵣ ( ̃a_{i+1,j}/ h_{i+1} - (h_i+h_{i+1})/(h_i*h_{i+1}) ̃a_{ij} +  a_{i-1,j}/ h_{i})

    A = zeros(BigFloat, 2N - 1, 2N - 1)

    Threads.@threads for i = 1:2N-1
        A[i, :] = 2 * CR * (A_1[i+1, :] / h[i+1] / (h[i] + h[i+1]) - A_1[i, :] / (h[i] * h[i+1]) + A_1[i-1, :] / h[i] / (h[i] + h[i+1]))
    end




    # truncation error

    U = u_exact.(x, α)
    F = parent(A) * U[1:2N-1]

    # plot(x[1:2N-1], F)

    U_s = parent(A) \ ones(2N - 1)
    # plot(x[1:2N-1], U[1:2N-1]-U_s)
    # plot!(x[1:2N-1], U_s)

    Si = sum(parent(A), dims=2)

    return x, U, F, U_s, Si
end

Ns = [100, 200, 400, 800]
num = length(Ns)
# trunc_err = zeros(num)
u_err = zeros(num)
RSi = zeros(num)
RS1 = zeros(num)
RSN = zeros(num)


# S = zeros(num)



Threads.@threads for l = 1:num
    @show N = Int64(Ns[l])
    (x, U, F, U_s, Si) = num_err(Ω, r, α, N)
    # trunc_err[l] = abs.(F .- 1) |> parent |> maximum
    u_err[l] = abs.(U[1:2N-1] - U_s) |> maximum
    # R = parent(abs.(F .- 1)) .* (x[1:2N-1] .^ α)
    R = parent(abs.(F .- 1)) ./ Si

    # println(Ri[l])
    RSi[l] = maximum(R[1:N])
    RS1[l] = (F[1] - 1) / Si[1]
    RSN[l] = (F[N] - 1) #/ Si[N]

    # S[l] = maximum( abs.(Si[1:N] .* x[1:N] .* α) )
    println("N = ", Ns[l], "  over")
end

# trunc_err, u_err

# plot(trunc_err)
# plot!(u_err)
# plot(u_err .* Ns.^(α/2))
# plot(Ri)
# plot(Ri .* Ns .^ (α / 2))
# plot(S)

@show diff( abs.(RSi) .|> log)./diff(log.(Ns))
@show diff( abs.(RS1) .|> log)./diff(log.(Ns))
@show diff( abs.(RSN) .|> log)./diff(log.(Ns))

@show diff( abs.(u_err) .|> log)./diff(log.(Ns))

u_err
