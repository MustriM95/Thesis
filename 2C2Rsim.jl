cd("C:\\Users\\micho\\github\\Thesis")

using DifferentialEquations
using LinearAlgebra
using Plots
using ModelingToolkit
using Distributions
using Random

## We define the dimension of our system, number of consumers and number
## of resources

M = 2
N = 3

## We initialize the parameters: uptake matrix, leakage matrix, inflow, and
## maitenance biomass

u = zeros(N, M)
l = zeros(N, M, M)
ρ = zeros(M, 1)
m = zeros(N, 1)



## Define preliminary variables to build uptake matrix and populate uᵢα
θ = zeros(N, M)
Ω = zeros(N, 1)
T = zeros(N, 1)
L = 0.4

## Construct uniform and beta distribution to draw parameters from
dU = Uniform(0,1)
dB = Beta(10, 20)

## Define a kroenecker delta function
δ(x, y) = ==(x, y)

## Sample concentration parameters for each consumer from uniform distribution
for i in 1:N
    θ[i, :] = rand(dU, M)
end

## Sample specialisation parameter for each consumer
Ω = 100 * rand(dU, N)

## Sample total uptake capacity per consumer
T = rand(dB, N)

## Generate uptake matrix from a dirichlet distribution
for i in 1:N
    dD = Dirichlet(Ω[i]*θ[i,:])
    u[i, :] = rand(dD) * T[i]
end

## Generate leakage tensor from dirichlet distribution
for i in 1:N
    dD = Dirichlet(Ω[i]*θ[i,:])
    for α in 1:M
        l[i, α, :] = rand(dD) * L
    end
end

## Sample inflow parameter from uniform distribution
ρ = rand(M)

## sample maitenance parameter from uniform distribution
m = rand(N)

## Generate ODE equations
@named sys1 = MiCRM(u = u, l = l, m = m, ρ = ρ, M = M, N = N)

## Initial conditions matrix
u0 = fill(1, N+M)
u0 = [states(sys1)[i] => u0[i] for i = eachindex(u0)]

## Generate ODE problem
prob = ODEProblem(sys1, u0, (0.0, 100.0), [], jac = true)

## Solve ODEs
sol = solve(prob)

## Plot solutions
plot(sol)

###############################################################################
# Functions
###############################################################################

function MiCRM(; name, u, l, ρ, m, M, N)
    M = M
    N = N
    ps = @parameters l[1:N, 1:M, 1:M] = l u[1:N, 1:M] = u ρ[1:M] = ρ m[1:N] = m
    sts = @variables t C[1:N](t) R[1:M](t)
    D = Differential(t)
    l = collect(l)
    u = collect(u)
    ρ = collect(ρ)
    m = collect(m)
    C = collect(C)
    R = collect(R)

    eqns = []

    for i in 1:N
        RHS = -C[i] * m[i]
        for α in 1:M
            RHS += C[i] * R[α] * u[i, α]
            for β in 1:M
                RHS += -C[i] * R[α] * u[i, α] * l[i, α, β]
            end
        end
        push!(eqns, D(C[i]) ~ RHS)
    end

    # Build system of equations for resources

    for β in 1:M
        RHS = ρ[β]
        for i in 1:N
            RHS += -u[i, β] * R[β] * C[i]
            for α in 1:M
                RHS += l[i, α, β] * u[i, α] * C[i] * R[α]
            end
        end
        push!(eqns, D(R[β]) ~ RHS)
    end

    ODESystem(eqns, t; name)
end

function Eff_LV(; name, u, l, ρ, m, M, N, Ceq, Req)
    A = zeros(M, M)
    K = zeros(N)
    r = zeros(N)
    ℵ = zeros(N, N)
    for α in 1:M
        for β in 1:M
            for i in 1:N
                A[α, β] += l[i, α, β] * u[i, β] * Ceq[i] -
                u[i, β] * Ceq[i] * δ(α, β)
            end
        end
    end


end

sol[1, 35]
sol[2, 35]
sol[3, 35]
sol[4, 35]
sol[5, 35]

Ceq = sol[1:N, 35]
A = zeros(M, M)

for α in 1:M
    for β in 1:M
        for i in 1:N
            A[α, β] += l[i, α, β] * u[i, β] * Ceq[i] -
            u[i, β] * Ceq[i] * δ(α, β)
        end
    end
end
