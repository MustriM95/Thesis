cd("C:\\Users\\micho\\github\\Thesis")

using DifferentialEquations
using LinearAlgebra
using Plots
using ModelingToolkit
using Distributions
using Random
using ForwardDiff
using DiffEqCallbacks
using StaticArrays


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


M = 2
N = 2

noise = Normal(0.001, 0.01)
NoiseM = rand(noise, N, M)

θ = 1.0*I(N)

θ = θ + broadcast(abs, NoiseM)

##@named p = MiCRM_par_therm(N = N, M = M, d = 0.05, μ = 0.2, η = 0.25, θ = θ)
@named p = MiCRM_par(N = N, M = M, L=0.7, μ = 0.2, θ = θ)

@named sys1 = MiCRM(p = p)

@named sys1 = update_parameters!(sys1, p_new = p)

u0 = fill(states(sys1)[1] => 0.0, (N+M))
for i in 1:N
    u0[i] = states(sys1)[i] => 0.1
end

for α in (N+1):(N+M)
    u0[α] = states(sys1)[α] => 0.1
end

prob = ODEProblem(sys1, u0, (0.0, 1000.0), [])

cb = AutoAbstol()

sol =solve(prob, TRBDF2(),reltol=1e-9, abstol=1e-9, saveat=1, callback = cb)
plot(sol)


@named Jac = MiCRM_jac(p=p, symbolic = false, sol=sol)

tc = zeros(N)
tR = zeros(M)

[tc[i]  = 1/abs(Jac[i, i]) for i in 1:N]
[tR[α] = 1/abs(Jac[α+N, α+N]) for α in 1:M]

trc = zeros(N, M)
for i in 1:N
    for α in 1:M
        trc[i, α] = tR[α]/(tc[i] + tR[α])
    end
end

trc_min = minimum(filter(!isnan,trc))

@named LV1 = Eff_LV_params(p = p, sol = sol, verbose = true)

LV1[:A]
LV1[:∂R]
LV1[:ℵ]
LV1[:r]

@named LVM = EFF_LV_sys(p_lv = LV1)

u0 = fill(0.1, N)
u0 = [states(LVM)[i] => u0[i] for i = eachindex(u0)]

prob_LV = ODEProblem(LVM, u0, (0.0, 500.0))

sol_LV = solve(prob_LV, TRBDF2(), reltol=1e-9, abstol=1e-9, saveat=1, callback = cb)

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

save("MCRM_fail11.jld", "data", p)


p = load("MCRM_fail5.jld", "data")
