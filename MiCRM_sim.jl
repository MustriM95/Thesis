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
N = 2

include("MiCRM_par.jl")
include("MiCRM_ODESys.jl")
include("Eff_Lv_params_alt.jl")
include("Eff_Lv_sys_alt.jl")

export MiCRM, MiCRM_par, Eff_Lv_params, Eff_Lv_sys

@named p = MiCRM_par(N = N, M = M)

@named sys1 = MiCRM(p = p)

u0 = fill(1, N+M)
u0 = [states(sys1)[i] => u0[i] for i = eachindex(u0)]

prob = ODEProblem(sys1, u0, (0.0, 100.0), [], jac = true)

sol = solve(prob)

plot(sol)

@named LV1 = Eff_LV_params(p = p, sol = sol)

@named LVM = EFF_LV_sys(p_lv = LV1)

u0 = fill(1, N)
u0 = [states(sys1)[i] => u0[i] for i = eachindex(u0)]

prob_LV = ODEProblem(LVM, u0, (0.0, 100.0), jac = true)

sol_LV = solve(prob_LV, Rosenbrock23())

plot(sol_LV)


p[:u]
p[:ω]
