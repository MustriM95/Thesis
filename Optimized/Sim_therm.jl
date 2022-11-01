("/home/michael/github*Thesis")

using Distributions
using LinearAlgebra
using DifferentialEquations
using Plots
using Sundials
using BenchmarkTools


include("Eff_LV_p_opt.jl")
include("LV_dx.jl")
include("Eff_LV_jac_opt.jl")
include("eq_t_cal.jl")
include("micrms_params_therm_opt.jl")
include("dx_therm.jl")

function F_m(N, M, kw)
    m = fill(0.2, N)
    return m
end

function F_ρ(N, M, kw)
    ρ = fill(0.2, M)
    return ρ
end

function F_ω(N, M, kw)
    ω = fill(0.0, M)
    return ω
end

N=5
M=5

p = MiCRM_par_therm(N, M; f_m=F_m, f_ρ=F_ρ, f_ω=F_ω)

x0 = fill(0.0, (N+M))
for i in 1:N
    x0[i] = 0.1
end

for α in (N+1):(N+M)
    x0[α] = 0.1
end

tspan = (0.0, 150.0)

## include("dx.jl")

prob = ODEProblem(dxx!, x0, tspan, p)

sol =solve(prob, CVODE_BDF(), saveat = 1)

plot(sol, idxs=[1, 2, 3, 4, 5])
plot(sol, idxs=[6, 7, 8, 9, 10])
