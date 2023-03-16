cd("/home/michael/github*Thesis")
cd("C:\\Users\\micho\\github\\Thesis")

using Distributions
using LinearAlgebra
using DifferentialEquations
using Plots
using Sundials
using BenchmarkTools

include("micrm_params.jl")
include("dx_v2.jl")
include("Eff_LV_p_opt.jl")
include("LV_dx.jl")
include("MiCRM_jac_opt.jl")
include("Eff_LV_jac_opt.jl")
include("MiCRM_test_opt.jl")
include("eq_t_cal.jl")
include("C:\\Users\\micho\\github\\Thesis\\Functions\\g_react.jl")


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

function F_u(N, M, kw)

    dU = Uniform(0, 1)

    u = zeros(N, M)

    if !haskey(kw, :θ)
        θ = zeros(N, M)

        for i in 1:N
            θ[i, :] = rand(dU, M)
        end
    else
        θ = kw[:θ]
    end


    if !haskey(kw, :Ω)
        ## Sample specialisation parameter for each consumer
        Ω = fill(1.0, N)
    else
        Ω = kw[:Ω]
    end

    ## Sample total uptake capacity per consumer
    T = fill(1.0, N)

    ## Generate uptake matrix from a dirichlet distribution
    for i = 1:N
        dD = Dirichlet(Ω[i] * θ[i, :])
        u[i, :] = rand(dD) * T[i]
    end
    return u

end

function F_l(N, M, kw)
    l = zeros(N, M, M)

    if !haskey(kw, :sim)
        sim = 1.0
    else
        sim = kw[:sim]
    end

    if !haskey(kw, :θ)
        dU = Uniform(0, 1)
        ϕ = rand(dU, M)
    else
        θ = kw[:θ]
    end

    θ_uni = zeros(N)
        for i in 1:N
            θ_uni += θ[i, :]/N
        end
    θ_uni
    min_index = 1
        
    for i in 2:N
        temp = θ_uni[1]
        if θ_uni[i] < temp
            min_index = i
        end
    end
        
    ortho = zeros(N)
    ortho[min_index] += θ_uni[min_index]

    ϕ = fill(0.0, M)
    ϕ = θ_uni*sim + (1/sim)*ortho
    ϕ = ϕ/norm(ϕ)

    dD = Dirichlet(100*ϕ[:])

    for i = 1:N
        for α = 1:M
            l[i, α, :] = rand(dD) * kw[:L]
        end
    end
    return l
end

EV = nothing

N=3
M=3

noise = Normal(1.01, 1.01)
NoiseM = rand(noise, N, M)
dU = Uniform(0, 1)

θ = 1.0*I(N)

θ = θ + broadcast(abs, NoiseM)

Ω = fill(100.0, N)

p = generate_params(N, M; f_u=F_u, f_m=F_m, f_ρ=F_ρ, f_ω=F_ω, L=0.3, θ=θ, Ω = Ω, sim=0.9)


x0 = fill(0.0, (N+M))
for i in 1:N
    x0[i] = 0.1
end

for α in (N+1):(N+M)
    x0[α] = 0.1
end

tspan = (0.0, 15000.0)

prob = ODEProblem(dxx!, x0, tspan, p)

sol =solve(prob, CVODE_BDF(), saveat = 1)

jac = MiCRM_jac(p=p, sol=sol)
if ==(EV, nothing)
    EV = complex(eigvals(jac))
else
    EV = append!(EV, eigvals(jac))
end

for α in (N+1):(N+M)
    x0[α] = 0.1
end

tspan = (0.0, 1500.0)

## include("dx.jl")

prob = ODEProblem(dxx!, x0, tspan, p)

sol =solve(prob, CVODE_BDF(), saveat = 1)

plot(sol, idxs=[1, 2, 3])
plot(sol, idxs=[4, 5, 6])

jac = MiCRM_jac(p=p, sol=sol)
g_react(J=jac)
#######################################################################################

p_lv = Eff_LV_params(p=p, sol=sol)

u0_LV = fill(0.0, N)
for i in 1:N
    u0_LV[i] = 0.1
end

LV_prob = ODEProblem(LV_dx!, u0_LV, tspan, p_lv)

sol_LV = solve(LV_prob, Tsit5(), reltol=1e-8, abstol=1e-8, saveat = 1)

plot(sol_LV)

LV_jac = Eff_Lv_Jac(p_lv=p_lv, sol=sol)

g_react(J=LV_jac)


eigen(LV_jac)

###############################################################################

eq_t = eq_t_cal(sol)

MAPE = 0.0
for i in 1:N
    for t in 1:Int(eq_t)
        MAPE += (100/(eq_t*N))*abs(sol[i, t] - sol_LV[i, t])/(sol[i, t])
    end
end

MPE = 0.0
for i in 1:N
    for t in 1:Int(eq_t)
        MPE += (100/(eq_t*N))*(sol[i, t] - sol_LV[i, t])/(sol[i, t])
    end
end

MAPE
MPE

SMAPE = 0.0
for i in 1:N
    for t in 1:Int(eq_t)
        SMAPE += log(10, abs(sol_LV[i, t]/sol[i, t]))/(eq_t*N)
    end
end

SMAPE

p_therm = MiCRM_par_therm(N, M)
x0 = fill(0.0, (N+M))
for i in 1:N
    x0[i] = 0.1
end

for α in (N+1):(N+M)
    x0[α] = 0.1
end

tspan = (0.0, 100.0)

## include("dx.jl")

prob = ODEProblem(dx_therm!, x0, tspan, p_therm)

sol =solve(prob, CVODE_BDF(), saveat = 1)

plot(sol, idxs=[1, 2, 3])
plot(sol, idxs=[4, 5, 6])