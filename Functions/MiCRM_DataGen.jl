cd("C:\\Users\\micho\\github\\Thesis")

using DifferentialEquations
using LinearAlgebra
using Plots
using ModelingToolkit
using Distributions
using Random
using ForwardDiff
using DiffEqCallbacks
using DataFrames

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

N = 5
M = 5

θ = zeros(N, M)

θ[1, :] = [1.0, 0.001, 0.001, 0.001, 0.001]
θ[2, :] = [0.001, 1.0, 0.001, 0.001, 0.001]
θ[3, :] = [0.001, 0.001, 1.0, 0.001, 0.001]
θ[4, :] = [0.001, 0.001, 0.001, 1.0, 0.001]
θ[5, :] = [0.001, 0.001, 0.001, 0.001, 1.0]

NoiseD = Normal(0.5, 0.2)
NoiseM = rand(NoiseD, N, M)

θ = θ + NoiseM

L = 0.3
μ = 0.2

@named p = MiCRM_par(N = N, M = M, L=L, μ = μ, θ = θ)

p[:u]

mean(cor(p[:u]) - Matrix(1.0I, N, M))

u_tr = zeros(N, M)

for i in 1:N
    for α in 1:M
        if i != α
            u_tr[i, α] = 0.0
        else
            u_tr[i, α] = p[:u][i, α]
        end
    end
end

off_diag_u = p[:u] - u_tr

sqrt(M*var(off_diag_u))

@named sys1 = MiCRM(p = p)

u0 = fill(states(sys1)[1] => 0.0, (N+M))
for i in 1:N
    u0[i] = states(sys1)[i] => 0.1
end

for α in (N+1):(N+M)
    u0[α] = states(sys1)[α] => 0.1
end


prob = ODEProblem(sys1, u0, (0.0, 700.0), [], jac = true)

sol =solve(prob,reltol=1e-8,abstol=1e-8, saveat=1)

plot(sol, vars = [1, 2, 3, 4, 5])

@named Jac = MiCRM_jac(p=p, symbolic = false, sol=sol)

@named LV1 = Eff_LV_params(p = p, sol = sol, verbose = true)

LV1[:ℵ]
LV1[:r]

@named LVM = EFF_LV_sys(p_lv = LV1)

u0 = fill(0.1, N)
u0 = [states(LVM)[i] => u0[i] for i = eachindex(u0)]

prob_LV = ODEProblem(LVM, u0, (0.0, 700.0), jac = true)

sol_LV = solve(prob_LV, reltol=1e-8,abstol=1e-8, saveat=1)

plot(sol_LV)

C_eq = zeros(N)
C_LV_eq = zeros(N)

for i in 1:N
    C_eq[i] = sol[i, 300]
end

for i in 1:N
    C_LV_eq[i] = sol_LV[i, 300]
end

@named LV_Jac = Eff_Lv_Jac(p_lv = LV1, sol = C_eq)

Eig = eigvals(Jac)

EigLV = eigvals(LV_Jac)

domEig = maximum(real(Eig))
domEigLV = maximum(real(EigLV))

diff_m_p = zeros(N)
diff_m_p = [minimum(broadcast(abs, sol[i, :] - sol_LV[i, :])) for i in 1:N]

null = fill(0.1, length(sol))
diff_m_n = zeros(N)
diff_m_n = [minimum(broadcast(abs, sol[i, :] - null)) for i in 1:N]

logL_p = 0.0

for j in 1:N
    av_sum = 0.0
    for i in 1:length(sol)
        av_sum += abs(sol[j, i] - sol_LV[j, i])/(sol[j, i])
    end
    logL_p = av_sum/N
end

logL_null = 0.0

for j in 1:N
    av_sum = 0.0
    for i in 1:length(sol)
        av_sum += abs(sol[j, i] - 0.1)/(sol[j, i])
    end
    logL_null += av_sum/N
end


p_R2 = 1 - (log(logL_p)/log(logL_null))

Eq_chi = 0.0

for i in 1:N
    Eq_chi += ((C_eq[i] - C_LV_eq[i])^2)
end

Eq_chi

u_tr = zeros(N, M)

for i in 1:N
    for α in 1:M
        if i != α
            u_tr[i, α] = 0.0
        else
            u_tr[i, α] = p[:u][i, α]
        end
    end
end

off_diag_u = p[:u] - u_tr

var(off_diag_u)
