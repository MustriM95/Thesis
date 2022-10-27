"""
This function creates a random set of parameters for the MiCRM.
Parameters are constrained by thermal dissipation d.

"""

function def_l(N, M, kw)
    L = zeros(N, M)
    l = zeros(N, M, M)

    for i in 1:N
        for α in 1:M
            L[i, α] = T[i]*(1-η[i, α]-d[i, α])
        end
    end

    ## Generate leakage tensor from dirichlet distribution
    for i in 1:N
        dD = Dirichlet(Ω[i]*θ[i,:])
        for α in 1:M
            l[i, α, :] = rand(dD)*L[i, α]
        end
    end


function MiCRM_par_therm(N, M; f_u=def_u, f_l=def_l, f_ρ=def_ρ, f_ω=def_ω, f_m=def_m, kwargs...)
    ## We initialize the parameters: uptake matrix, leakage tensor, inflow, and
    ## maitenance biomass

    kw = Dict{Symbol, Any}(kwargs)

    u = zeros(N, M)
    ρ = zeros(M, 1)
    ω = zeros(M, 1)
    m = zeros(N, 1)

    ## Define preliminary variables to build uptake matrix and populate uᵢα
    Ω = zeros(N, 1)
    T = zeros(N, 1)

    if d == nothing
        d = rand(N, M)
    else
        dNorm = Normal(d, d*0.1)
        d = rand(dNorm, N, M)
    end

    if η == nothing
        η = rand(N, M)
    else
        ηNorm = Normal(η, η*0.1)
        η = rand(ηNorm, N, M)
    end

    ## Construct uniform and beta distribution to draw parameters from
    dU = Uniform(0,1)
    dB = Beta(4, 3)
    dN = Normal(μ, 0.1*μ)
    dΩ = Uniform(1, 10000)


    ## Sample concentration parameters for each consumer from uniform distribution
    if θ == nothing
        θ = zeros(N, M)

        for i in 1:N
            θ[i, :] = rand(dU, M)
        end
    end

    ## Sample specialisation parameter for each consumer
    Ω = rand(dΩ, N)

    ## Sample total uptake capacity per consumer
    T = rand(dB, N)

    ## Generate uptake matrix from a dirichlet distribution
    for i in 1:N
        dD = Dirichlet(Ω[i]*θ[i, :])
        u[i, :] = rand(dD)*T[i]
    end

## Calculation of total row sum leakage from efficiencyt and dissipation
# terms

    for i in 1:N
        for α in 1:M
            L[i, α] = T[i]*(1-η[i, α]-d[i, α])
        end
    end




    ## Generate leakage tensor from dirichlet distribution
    for i in 1:N
        dD = Dirichlet(Ω[i]*θ[i,:])
        for α in 1:M
            l[i, α, :] = rand(dD)*L[i, α]
        end
    end

    ## Sample inflow parameter from uniform distribution
    ρ = rand(M)*(1/M)
    ω = rand(M)*(1/(2*M))

    ## sample maitenance parameter from uniform distribution
    m = [rand(dN)*T[i]^(2/3) for i in 1:N]

    p = Dict(:l => l, :ρ => ρ, :ω => ω, :m => m, :M => M, :N => N, :u => u)
end
