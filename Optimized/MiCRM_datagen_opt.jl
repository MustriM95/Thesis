cd("/home/michael/github*Thesis")

using DifferentialEquations
using LinearAlgebra
using Plots
using Distributions
using Random
using DiffEqCallbacks
using DataFrames
using CSV
using BenchmarkTools


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
        Ω = fill(kw[:Ω], N)
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

## We define the dimension of our system, number of consumers and number
## of resources

include("micrm_params.jl")
include("dx_v2.jl")
include("Eff_LV_p_opt.jl")
include("LV_dx.jl")
include("MiCRM_jac_opt.jl")
include("Eff_LV_jac_opt.jl")
include("MiCRM_test_opt_v2.jl")

t_span = 100000.0

N = 14
M = 14

θ = zeros(N, M)

Ω = fill(1.0, N)

itr = 0
results = nothing
Ω_vec = [1.0, 10, 100, 1000]
σ_vec = [0.1, 0.15, 0.2, 0.3, 0.4, 0.5]
μ_vec = [0.1, 0.15, 0.2, 0.3, 0.4, 0.5]

for μ in μ_vec
    for σ in σ_vec
        for Ω in Ω_vec
            L = 0.0
            while L < 0.9
                Sim = MiC_test(μ=μ, σ=σ, L=L,  N=N, M=M, Ω = Ω, t_span=t_span)
                if itr == 0
                    results= N, M, L, μ, σ, Sim[:NO], Sim[:CO], Sim[:SMAPE], Sim[:Eq_SMAPE], Sim[:ℵ_m], Sim[:r_m],
                     Sim[:eq_t], Sim[:domEig], Sim[:domEigLV], Sim[:C_sur], Sim[:trc_max]
                else
                    temp = N, M, L, μ, σ, Sim[:NO], Sim[:CO],  Sim[:SMAPE], Sim[:Eq_SMAPE], Sim[:ℵ_m], Sim[:r_m],
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

C_names = [:N, :M ,:Leakage, :Noise_m, :Noise_std, :NO, :CO, :SMAPE, :Eq_SMAPE, :I_m, :R_m,
 :eq_t, :domEig, :domEigLV, :C_sur, :trc_max]



df = DataFrame(results)
df = rename(df, C_names)

CSV.write("./Data/14x14_opt_v2.csv", df)
