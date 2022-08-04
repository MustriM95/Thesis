cd("/home/pawarlab/Thesis/HPC")
using Distributed

@everywhere using DifferentialEquations
@everywhere using LinearAlgebra
@everywhere using Plots
@everywhere using ModelingToolkit
@everywhere using Distributions
@everywhere using Random
@everywhere using ForwardDiff
@everywhere using DiffEqCallbacks
@everywhere using DataFrames
@everywhere using CSV
@everywhere using Distributed
@everywhere using SharedArrays



## We define the dimension of our system, number of consumers and number
## of resources

@everywhere include("MiCRM_par.jl")
@everywhere include("MiCRM_ODESys.jl")
@everywhere include("Eff_Lv_params_alt.jl")
@everywhere include("Eff_Lv_sys_alt.jl")
@everywhere include("MiCRM_jac.jl")
@everywhere include("Eff_Lv_jac.jl")
@everywhere include("MiCRM_par_therm.jl")
@everywhere include("MiC_test.jl")
@everywhere include("Sim_call.jl")


export MiCRM, MiCRM_par, Eff_Lv_params, Eff_Lv_sys, MiCRM_jac, Eff_Lv_Jac
export MiCRM_par_therm, MiC_test, Sim_call

tspan = 10000.0

N = 3
M = 3

θ = 1.0*I(N)

Ω = fill(1.0, N)

μ = 0.0
σ = 0.0
@everywhere itr = SharedVector{Int64}(1)
itr[1] = 1
results = SharedArray{Float64}(16, 1440)
Ω_vec = [1.0, 10, 100, 1000]
σ_vec = [0.001, 0.1, 0.2, 0.3, 0.4, 0.5]
μ_vec = [0.001, 0.1, 0.2, 0.3, 0.4, 0.5]
L_vec = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]


@sync @distributed for L in L_vec
    for σ in σ_vec
        for Ω in Ω_vec
            for μ in μ_vec
                Sim = MiC_test(μ=μ, σ=σ, L=L,  N=N, M=M, θ=θ, Ω = Ω, t_span=tspan)
                results[:, itr] = [convert(Float64, N), convert(Float64, M), L, μ, σ, Sim[:NO], Sim[:CO],
                 Sim[:MSE], Sim[:Eq_MSE], Sim[:ℵ_m], Sim[:r_m],convert(Float64, Sim[:eq_t]),
                 Sim[:domEig], Sim[:domEigLV], Sim[:C_sur], Sim[:trc_max]]
                print(itr)
                itr[1] += 1
            end
        end
    end
end


C_names = [:N, :M ,:Leakage, :Noise_m, :Noise_std, :NO, :CO, :MSE, :Eq_MSE, :I_m, :R_m,
 :eq_t, :domEig, :domEigLV, :C_sur, :trc_max]

res = convert(Array, results)
df = DataFrame(res', :auto)
df = rename(df, C_names)

CSV.write("parallel_text.csv", df)
