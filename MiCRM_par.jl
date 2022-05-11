"""
This function creates a random set of parameters for the MiCRM.
Alternatively, the option of structuring certain parameters will
be implemented.

"""

function MiCRM_par(; name, N, M)
    ## We initialize the parameters: uptake matrix, leakage matrix, inflow, and
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
    L = 0.1

    ## Construct uniform and beta distribution to draw parameters from
    dU = Uniform(0,1)
    dB = Beta(10, 20)

    ## Sample concentration parameters for each consumer from uniform distribution
    for i in 1:N
        θ[i, :] = rand(dU, M)
    end

    ## Sample specialisation parameter for each consumer
    Ω = 100 * rand(dU, N)

    ## Sample total uptake capacity per consumer
    T = rand(dB, N)

    ## Generate uptake matrix from a dirichlet distribution
    for i in 1:N
        dD = Dirichlet(Ω[i]*θ[i,:])
        u[i, :] = rand(dD) * T[i]
    end

    ## Generate leakage tensor from dirichlet distribution
    for i in 1:N
        dD = Dirichlet(Ω[i]*θ[i,:])
        for α in 1:M
            l[i, α, :] = rand(dD) * L
        end
    end

    ## Sample inflow parameter from uniform distribution
    ρ = rand(M)
    ω = rand(M)*0.5

    ## sample maitenance parameter from uniform distribution
    m = rand(N)

    p = Dict(:l => l, :ρ => ρ, :ω => ω, :m => m, :M => M, :N => N, :u => u)
end
