cd("C:\\Users\\micho\\github\\Thesis")

using DifferentialEquations
using LinearAlgebra
using Plots
using ModelingToolkit
using Distributions
using Random
using ForwardDiff
using DiffEqCallbacks

## We define the dimension of our system, number of consumers and number
## of resources

include("MiCRM_par.jl")
include("MiCRM_ODESys.jl")
include("Eff_Lv_params_alt.jl")
include("Eff_Lv_sys_alt.jl")
include("MiCRM_jac.jl")
include("Eff_Lv_jac.jl")
include("MiCRM_par_therm.jl")


export MiCRM, MiCRM_par, Eff_Lv_params, Eff_Lv_sys, MiCRM_jac, Eff_Lv_Jac
export MiCRM_par_therm


M = 3
N = 3

noise = Normal(0.1, 0.01)


θ = zeros(N, M)

θ[1, : ] = [1.0, 0.01, 0.01]
θ[2, : ] = [0.01, 1.0, 0.01]
θ[3, : ] = [0.01, 0.01, 1.0]
##θ[4, : ] = [1., 1., 1.]
##θ[5, : ] = [1., 0.0001, 1.]

##@named p = MiCRM_par_therm(N = N, M = M, d = 0.05, μ = 0.2, η = 0.25, θ = θ)
@named p = MiCRM_par(N = N, M = M, L=0.5, μ = 0.2, θ = θ)


cor(p[:u][1, :], p[:u][2, :])
cor(p[:u][1, :], p[:u][3, :])
cor(p[:u][2, :], p[:u][3, :])

p[:ρ] = [0.2, 0.2, 0.2]
p[:ω] = zeros(M)

@named sys1 = MiCRM(p = p)

u0 = fill(states(sys1)[1] => 0.0, (N+M))
for i in 1:N
    u0[i] = states(sys1)[i] => 0.1
end

for α in (N+1):(N+M)
    u0[α] = states(sys1)[α] => 0.1
end


prob = ODEProblem(sys1, u0, (0.0, 500.0), [], jac = true)

sol =solve(prob,reltol=1e-8,abstol=1e-8, saveat=1)
plot(sol)


@named Jac = MiCRM_jac(p=p, symbolic = false, sol=sol)

ct = [1/abs(Jac[i, i]) for i in 1:(N+M)]

Jᵣ = zeros((N+M), (N+M))

[Jᵣ[i, :] = Jac[i, :] for i in 1:(N+M)]


maximum([1/abs(Jᵣ[1, i]) for i in 1:(N+M)])


@named LV1 = Eff_LV_params(p = p, sol = sol, verbose = true)

LV1[:A]
LV1[:∂R]
LV1[:ℵ]
LV1[:r]

@named LVM = EFF_LV_sys(p_lv = LV1)

u0 = fill(0.1, N)
u0 = [states(LVM)[i] => u0[i] for i = eachindex(u0)]

prob_LV = ODEProblem(LVM, u0, (0.0, 500.0), jac = true)

sol_LV = solve(prob_LV)

plot(sol_LV)

@named LV_Jac = Eff_Lv_Jac(p_lv = LV1, sol = sol_LV)

ct_LV = [1/abs(LV_Jac[i, i]) for i in 1:(N)]

eigvals(Jac)

eigvals(LV_Jac)

[LV_Jac[i, i] for i in 1:N]

Ceq = sol[1:N, length(sol)]

Ceq_LV = sol_LV[1:N, length(sol_LV)]

Ceq - Ceq_LV

## Saving parameter dictionary
using JLD, HDF5

save("MCRM_fail10.jld", "data", p)


p = load("MCRM_fail5.jld", "data")
