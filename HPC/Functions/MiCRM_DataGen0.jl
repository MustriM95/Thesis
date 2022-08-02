cd("C:\\Users\\micho\\github\\Thesis\\HPC")

using DifferentialEquations
using LinearAlgebra
using Plots
using ModelingToolkit
using Distributions
using Random
using ForwardDiff
using DiffEqCallbacks
using DataFrames
using CSV


## We define the dimension of our system, number of consumers and number
## of resources

include("MiCRM_par.jl")
include("MiCRM_ODESys.jl")
include("Eff_Lv_params_alt.jl")
include("Eff_Lv_sys_alt.jl")
include("MiCRM_jac.jl")
include("Eff_Lv_jac.jl")
include("MiCRM_par_therm.jl")
include("MiC_test.jl")


export MiCRM, MiCRM_par, Eff_Lv_params, Eff_Lv_sys, MiCRM_jac, Eff_Lv_Jac
export MiCRM_par_therm, MiC_test

t_span = 10000.0

N = 2
M = 2

θ = 1.0*I(N)

Ω = fill(1.0, N)

μ = 0.0
σ = 0.0

global itr = 0
global results = nothing
Ω_vec = [1.0, 10, 100, 1000]
σ_vec = [0.001, 0.1, 0.2, 0.3, 0.4, 0.5]
μ_vec = [0.001, 0.1, 0.2, 0.3, 0.4, 0.5]

for μ in μ_vec
    for σ in σ_vec
        for Ω in Ω_vec
            L = 0.0
            while L < 0.9
                Sim = MiC_test(μ=μ, σ=σ, L=L,  N=N, M=M, θ=θ, Ω = Ω, t_span=t_span)
                if itr == 0
                    results= N, M, L, μ, σ, Sim[:NO], Sim[:MSE], Sim[:Eq_MSE], Sim[:ℵ_m], Sim[:r_m],
                     Sim[:eq_t], Sim[:domEig], Sim[:domEigLV], Sim[:C_sur], Sim[:trc_max]
                else
                    temp = N, M, L, μ, σ, Sim[:NO], Sim[:MSE], Sim[:Eq_MSE], Sim[:ℵ_m], Sim[:r_m],
                     Sim[:eq_t], Sim[:domEig], Sim[:domEigLV], Sim[:C_sur], Sim[:trc_max]
                    results = vcat(results, temp)
                end
                itr += 1
                L += 0.1
                println(itr)
            end
        end
    end
end

C_names = [:N, :M ,:Leakage, :Noise_m, :Noise_std, :NO, :MSE, :Eq_MSE, :I_m, :R_m,
 :eq_t, :domEig, :domEigLV, :C_sur, :trc_max]

df = DataFrame(results)
df = rename(df, C_names)

CSV.write("16x16.csv", df)
