cd("C:\\Users\\micho\\github\\Thesis")

using Distributions, LinearAlgebra, DataFrames, CSV
using DifferentialEquations, DiffEqCallbacks
include("micrm_params.jl")
include("par_functions.jl")
include("MiCRM_test_opt_v3.jl")
include("dx_v2.jl")
include("Eff_LV_p_opt.jl")
include("LV_dx.jl")
include("MiCRM_jac_opt.jl")
include("Eff_LV_jac_opt.jl")
include("C:\\Users\\micho\\github\\Thesis\\Functions\\g_react.jl")

function NiOv(; p)
    cossim(x, y) = dot(x, y) / (norm(x)*norm(y))

    NO = 0.0
    for i in 1:N
        for j in (i+1):N
            NO += cossim(p.u[i, :], p.u[j, :])*(2/(N*(N-1)))
        end
    end

    return NO
end

function CoOp(; p)
    eff_L = zeros(N, M)

    for j in 1:N
        for i in 1:M
            eff_L[j, :] += p.u[j, i]*p.l[j, i, :]
        end
    end

    cossim(x, y) = dot(x, y) / (norm(x)*norm(y))

    CO = 0.0
    for i in 1:N
        for j in (i+1):N
            CO += cossim(p.u[i, :], eff_L[j, :])*(2/(N*(N-1)))
        end
    end

    return CO
end

sim_draw = [0.5, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6]

Nsim_draw = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]

L_draw = [0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8]

N=4
M=4

noise = Normal(1.0, 1.0)
NoiseM = rand(noise, N, M)  

itr = 0
community = nothing
for L in L_draw
    for N_sim in Nsim_draw
        for sim in sim_draw

            θ = 1.0*I(N)

            θ_het = fill(1.0, N, M)

            θ = 100*N_sim*θ + broadcast(abs, NoiseM) + 30*(1.0 - N_sim)*θ_het

            Ω_mod = 10*(1/N_sim) + 100*N_sim

            Ω = fill(Ω_mod, N)

            p = generate_params(N, M; f_u=F_u, f_m=F_m, f_ρ=F_ρ, f_ω=F_ω, f_l=F_l, L=L, θ=θ, Ω = Ω, sim = sim)

            NO = NiOv(p=p)
            CO = CoOp(p=p)

            if itr == 0
                community = N, M, sim, N_sim, NO, CO, L
            else
                temp = N, M, sim, N_sim, NO, CO, L
                community = vcat(community, temp)
            end

            itr += 1

        end
    end
end

C_names = [:N, :M , :sim, :N_sim, :NO, :CO, :L]

df = DataFrame(community)
df = rename(df, C_names)
CSV.write("similarity_test.csv", df)

itr = 0
community_df = nothing
for L in L_draw
    for N_sim in Nsim_draw
        for sim in sim_draw

            θ = 1.0*I(N)

            θ_het = fill(1.0, N, M)

            θ = 100*N_sim*θ + broadcast(abs, NoiseM) + 30*(1.0 - N_sim)*θ_het

            Ω_mod = 10*(1/N_sim) + 100*N_sim

            Ω = fill(Ω_mod, N)

            p = generate_params(N, M; f_u=F_u, f_m=F_m, f_ρ=F_ρ, f_ω=F_ω, f_l=F_l, L=L, θ=θ, Ω = Ω, sim = sim)

            if itr == 0
                community_df = DataFrame([p])
            else
                temp = DataFrame([p])
                community_df = vcat(community_df, temp)
            end

            itr += 1

        end
    end
end

community_df

community_df

p_1 = NamedTuple(community_df[250, :])

MiC_test(p=p_1)

itr = 0
results = nothing
for row in eachrow(community_df)
    p = NamedTuple(row)
    Sim = MiC_test(p=p)

    if itr == 0
        results = p.N, p.M, p.L, Sim[:NO], Sim[:CO], Sim[:SMAPE], Sim[:Eq_SMAPE], 
         Sim[:eq_t], Sim[:domEig], Sim[:domEigLV], Sim[:C_sur], Sim[:trc_max], Sim[:gR], Sim[:gR_LV]
    else 
        temp = p.N, p.M, p.L, Sim[:NO], Sim[:CO], Sim[:SMAPE], Sim[:Eq_SMAPE], 
         Sim[:eq_t], Sim[:domEig], Sim[:domEigLV], Sim[:C_sur], Sim[:trc_max], Sim[:gR], Sim[:gR_LV]
        results = vcat(results, temp)
    end
    itr += 1

end

C_names = [:N, :M ,:Leakage, :NO, :CO, :SMAPE, :Eq_SMAPE,
 :eq_t, :domEig, :domEigLV, :C_sur, :trc_max, :gR, :gR_LV]



df = DataFrame(results)
df = rename(df, C_names)

CSV.write("4x4_newTest.csv", df)