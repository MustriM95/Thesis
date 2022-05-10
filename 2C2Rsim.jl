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
ω = zeros(M, 1)
m = zeros(N, 1)



## Define preliminary variables to build uptake matrix and populate uᵢα
θ = zeros(N, M)
Ω = zeros(N, 1)
T = zeros(N, 1)
L = 0.2

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
ω = rand(M)

## sample maitenance parameter from uniform distribution
m = rand(N)

## Generate ODE equations
@named sys1 = MiCRM(u = u, l = l, m = m, ρ = ρ, M = M, N = N, ω = ω)

## Initial conditions matrix
u0 = fill(1, N+M)
u0 = [states(sys1)[i] => u0[i] for i = eachindex(u0)]

## Generate ODE problem
prob = ODEProblem(sys1, u0, (0.0, 100.0), [], jac = true)

## Solve ODEs
sol = solve(prob)

## Plot solutions
plot(sol)

p = Dict(:u => u, :l => l, :ρ => ρ, :ω => ω, :m => m, :M => M, :N => N)

@named LV1 = Eff_LV_params(p = p, sol = sol)

@named LVM = EFF_LV_sys(p_lv = LV1)

u0 = fill(1, N)
u0 = [states(sys1)[i] => u0[i] for i = eachindex(u0)]

prob_LV = ODEProblem(LVM, u0, (0.0, 100.0), jac = true)

sol_LV = solve(prob_LV, Rosenbrock23())

plot(sol_LV)

Ceq = sol[1:N, length(sol)]
Ceq_LV = sol_LV[1:N, length(sol_LV)]

equations(LVM)

sol_LV[2, 100]

do something

###############################################################################
# Functions
###############################################################################

function MiCRM(; name, u, l, ρ, m, M, N, ω)
    M = M
    N = N
    ps = @parameters l[1:N, 1:M, 1:M] = l u[1:N, 1:M] = u ρ[1:M] = ρ m[1:N] = m
    sts = @variables t C[1:N](t) R[1:M](t)
    D = Differential(t)
    l = collect(l)
    u = collect(u)
    ρ = collect(ρ)
    ω = collect(ω)
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
        RHS = ρ[β] - R[β]*ω[β]
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



## Note: If total leakage is the same for all species and
## resources, the diagonal terms of the interaction matrix
## will be 1

function Eff_LV_params(; name, p, sol)
    @assert all([:l, :ρ, :ω, :m, :M, :N, :u] .∈ Ref(collect(keys(p)))) "missing parameters"

    M = p[:M]
    N = p[:N]
    l = p[:l]
    ρ = p[:ρ]
    ω = p[:ω]
    m = p[:m]
    u = p[:u]


    Ceq = sol[1:N, length(sol)]
    Req = sol[(N+1):(N+M), length(sol)]
    A = zeros(M, M)
    ∂R = zeros(M, N)
    ℵ = zeros(N, N)
    λ = zeros(N, M)
    O = zeros(N) # Dummy variable for first term of sum
    P = zeros(N) # Dummy variable for interaction term
    r = zeros(N) # Effective growth rates

## Calculating the A matrix from parameters at equilibrium
## solutions
    for α in 1:M
        for β in 1:M
            A[α, β] = -ω[α]
            for i in 1:N
                A[α, β] += l[i, α, β] * u[i, β] * Ceq[i] -
                u[i, β] * Ceq[i] * δ(α, β)
            end
        end
    end

## Calculating partial derivatives of eq resources with respect
## to consumers
    for α in 1:M
        for j in 1:N
            for β in 1:M
                for γ in 1:M
                    ∂R[α, j] += inv(A)[α, β] * u[j, β] * Req[γ]*(δ(β, γ)
                    - l[j, β, γ])
                end
            end
        end
    end

## Calculating total leakage of consumers per resource
    for i in 1:N
        for α in 1:M
            λ[i, α] = sum(l[i, α, :])
        end
    end

## Calculating interaction matrix
    for i in 1:N
        for j in 1:N
            for α in 1:M
                ℵ[i, j] += u[i, α]*(1 - λ[i, α]) * ∂R[α, j]
            end
        end
    end

## Calculating effective carrying capacities

    for i in 1:N
        for α in 1:M
            O[i] += u[i, α]*(1 - λ[i, α])*Req[α]
        end
    end

    for i in 1:N
        for j in 1:N
            P[i] += ℵ[i, j]*Ceq[j]
        end
    end

## Calculating effective intrinsic growth rates
    for i in 1:N
        r[i] = O[i] - P[i] -m[i]
    end

    return Dict(:ℵ => ℵ, :r => r, :N => N)
end

function EFF_LV_sys(; name, p_lv)
    @assert all([:ℵ, :r, :N] .∈ Ref(collect(keys(p_lv))))

    N = p_lv[:N]
    ℵ = ModelingToolkit.toparam(Symbolics.variables(:ℵ, 1:N, 1:N))
    r = ModelingToolkit.toparam(Symbolics.variables(:r, 1:N))

    sts = @variables t C[1:N](t)

    D = Differential(t)

    eqns = []

    for i in 1:N
        RHS = r[i]*C[i]
        for j in 1:N
            RHS += C[i]*ℵ[i, j]*C[j]
        end
        push!(eqns, D(C[i]) ~ RHS)
    end

    vars = vcat(ℵ[:], r)
    vals = vcat(p_lv[:ℵ][:], p_lv[:r])

    sys = ODESystem(eqns, t; name, defaults = Dict(vars[i] => vals[i] for i = eachindex(vars)))

    return(sys)
end
