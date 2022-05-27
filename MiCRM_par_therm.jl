"""
This function creates a random set of parameters for the MiCRM.
Parameters are constrained by thermal dissipation d.

"""

function MiCRM_par_therm(; name, N, M, η = nothing, d = nothing)
    ## We initialize the parameters: uptake matrix, leakage tensor, inflow, and
    ## maitenance biomass

    u = zeros(N, M)
    l = zeros(N, M, M)
    ρ = zeros(M, 1)
    ω = zeros(M, 1)
    m = zeros(N, 1)

    ## Define preliminary variables to build uptake matrix and populate uᵢα
    θ = zeros(N, M)
    Ω = zeros(N, 1)
    T = zeros(N, 1)
    L = zeros(N, M)

    d = rand(N, M)

    η = rand(N, M)


    ## Construct uniform and beta distribution to draw parameters from
    dU = Uniform(0,1)
    dB = Beta(4, 3)
    dN = Normal(0.5, 0.05)


    ## Sample concentration parameters for each consumer from uniform distribution
    for i in 1:N
        θ[i, :] = rand(dU, M)
    end

    ## Sample specialisation parameter for each consumer
    Ω = rand(dU, N)

    ## Sample total uptake capacity per consumer
    T = rand(dB, N)

    ## Generate uptake matrix from a dirichlet distribution
    for i in 1:N
        ϕ = zeros(M)
        ϕ = [1- θ[i, j] for j in 1:M]
        dD = Dirichlet(Ω[i]*ϕ[:])
        u[i, :] = rand(dD)*T[i]
    end

## Calculation of total row sum leakage from efficiencyt and dissipation
# terms

    for i in 1:N
        for α in 1:M
            L[i, α] = T[i]*(1-η[i, α])*(1-d[i, α])
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
    m = [rand(dN) for i in 1:N]

    p = Dict(:l => l, :ρ => ρ, :ω => ω, :m => m, :M => M, :N => N, :u => u)
end
