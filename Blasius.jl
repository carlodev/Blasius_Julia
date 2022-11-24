using Plots
using DifferentialEquations


function balsius_fun!(du, u, p, η)
    f, f1, f2 = u
    du[1] = f1
    du[2] = f2
    du[3] = -0.5*f2*f
end

function bc_blasius!(residual, u, p, η) # u[1] is the beginning of the time span, and u[end] is the ending
    residual[1] = u[1][1]   # the solution at the beginning of the time span should be -pi/2
    residual[2] = u[1][2]  # the solution at the end of the time span should be pi/2
    residual[3] = u[end][2] - 1.0
end
ηspan = (0,8)

blasius_problem = TwoPointBVProblem(balsius_fun!, bc_blasius!, [0.0, 0.0, 0.3321], ηspan)
balsius_solution = solve(blasius_problem, Shooting(Vern7()), dt=0.001) #dt in reality dη

plot(balsius_solution,  label = [ "f(η)" "f'(η)" "f''(η)"], legend= :topleft)
plot!(xlabel ="η")
savefig("Blasius.pdf")

f(η) = balsius_solution(η)[1]
fp(η) = balsius_solution(η)[2]
η(x,y) = y * (U₀ / (ν*x))^0.5
ψ(η) = (ν* U₀*x)*f(η)
u(x,y) = U₀ * fp(η(x,y))
v(x,y) =0.5 * (ν* U₀/x)^0.5 *(η(x,y) * fp(η(x,y))- f(η(x,y)))


# Testing
Reδ₀ = 795
ν = 1e-5
U₀ = 1.5 #freestream velocity
δ₀ = Reδ₀ * ν/ U₀ #The height of the boundary layer when it has reached 0.99U₀ 

#Position form the leading edge where to compute the physical boundary layer
x_in = (δ₀/5.29)^2 * U₀ / ν

yv = LinRange(0,δ₀,100)
u.(x_in, yv)

plot(yv, u.(x_in, yv))
plot(yv, v.(x_in, yv))
