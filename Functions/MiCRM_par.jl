"""
This function creates a random set of parameters for the MiCRM.
Alternatively, the option of structuring certain parameters will
be implemented.

"""

function MiCRM_par(; name, N, M, L, μ = 0.5, θ = nothing)
    ## We initialize the parameters: uptake matrix, leakage matrix, inflow, and
    ## maitenance biomass

    u = zeros(N, M)
    l = zeros(N, M, M)
    ρ = zeros(M, 1)
    ω = zeros(M, 1)
    m = zeros(N, 1)

    ## Define preliminary variables to build uptake matrix and populate uᵢα
    Ω = zeros(N, 1)
    T = zeros(N, 1)

    ## Construct uniform and beta distribution to draw parameters from
    dU = Uniform(0, 1)
    dB = Beta(4, 3)
    dN = Normal(μ, 0.1*μ)

    ## Sample concentration parameters for each consumer from uniform distribution
    if θ == nothing
        θ = zeros(N, M)

        for i in 1:N
            θ[i, :] = rand(dU, M)
        end
    end

    ## Sample specialisation parameter for each consumer
    Ω = fill(1, N)

    ## Sample total uptake capacity per consumer
    T = fill(1, N)

    ## Generate uptake matrix from a dirichlet distribution
    for i = 1:N
        dD = Dirichlet(Ω[i] * θ[i, :])
        u[i, :] = rand(dD) * T[i]
    end

    ## Generate leakage tensor from dirichlet distribution
    for i = 1:N
        dD = Dirichlet(Ω[i] * θ[i, :])
        for α = 1:M
            l[i, α, :] = rand(dD) * L
        end
    end

    ## Sample inflow parameter from uniform distribution
    ρ = fill(0.2, M)
    ω = zeros(M)

    ## sample maitenance parameter from normal distribution
    m = rand(dN, N)

    p = Dict(:l => l, :ρ => ρ, :ω => ω, :m => m, :M => M, :N => N, :u => u)
end
