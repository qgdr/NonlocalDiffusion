using Plots

f(x, y, α) = x*y *(1+ x^(2-α) + y^(2-α) ) + (1+x+y)^(3-α) + (1+x)*(1+y) * (1- (1+x)^(2-α) - (1+y)^(2-α) )

# n=101
# x = range(0, 2//10, length=n)
# y = x

# z = zeros(n,n)
# for i in 1:n, j in 1:n
#     z[i,j] = f(x[i], y[j], 1.1)
# end

# Plots.surface(x, y, z)

df_dx(x, y, α) = 1 + 2y + y^(3-α) - (1+y)^(3-α)  +  (3-α)* ( x^(2-α) * y + (1+x+y)^(2-α) -(1+x)^(2-α)*(1+y) )