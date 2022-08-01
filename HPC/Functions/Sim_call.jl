## Function for looping through simulations and extracting data points.

function Sim_call(M, N, tspan, θ)
    itr = 0
    results = nothing
    Ω_vec = [1.0, 10, 100, 1000]
    σ_vec = [0.001, 0.1, 0.2, 0.3, 0.4, 0.5]
    μ_vec = [0.001, 0.1, 0.2, 0.3, 0.4, 0.5]
    for μ in μ_vec
        for σ in σ_vec
            for Ω in Ω_vec
                L = 0.0
                while L < 0.9
                    Sim = @spawn MiC_test(μ=μ, σ=σ, L=L,  N=N, M=M, θ=θ, Ω = Ω, t_span=t_span)
                    wait(Sim)

                    if itr == 0
                        global results= N, M, L, μ, σ, Sim[:NO], Sim[:MSE], Sim[:Eq_MSE], Sim[:ℵ_m], Sim[:r_m],
                         Sim[:eq_t], Sim[:domEig], Sim[:domEigLV], Sim[:C_sur], Sim[:trc_max]
                    else
                        temp = N, M, L, μ, σ, Sim[:NO], Sim[:MSE], Sim[:Eq_MSE], Sim[:ℵ_m], Sim[:r_m],
                         Sim[:eq_t], Sim[:domEig], Sim[:domEigLV], Sim[:C_sur], Sim[:trc_max]
                        global results = vcat(results, temp)
                    end
                    global itr += 1
                    L += 0.1
                    println(itr)
                end
            end
        end
    end
    return results
end
