"""
This function creates a random set of parameters for the MiCRM.
Parameters are constrained by thermal dissipation d.

"""
def_ρ(N, M, kw) = ones(M)
def_ω(N, M, kw) = ones(M)

function def_m(N, M, kw)

    if !haskey(kw, :T)
        T = fill(1.0, N)
    else
        T = kw[:T]
    end

    m = [T[i]^(2/3) for i in 1:N]
    return m
end


function def_l(N, M, d, kw)
    L = zeros(N, M)
    l = zeros(N, M, M)
    effD = Uniform(0.1, 0.85)

    if !haskey(kw, :η)
        η = rand(effD, N, M)
    else
        η = kw[:η]
    end

    if !haskey(kw, :T)
        T = fill(1.0, N)
    else
        T = kw[:T]
    end

    for i in 1:N
        for α in 1:M
            L[i, α] = T[i]*(1-η[i, α]-d[i, α])
        end
    end

    ## Generate leakage tensor from dirichlet distribution
    ϕ = fill(1.0, M)
    dD = Dirichlet(ϕ[:])
    for i in 1:N
        for α in 1:M
            l[i, α, :] = rand(dD)*L[i, α]
        end
    end
    return l
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
        Ω = kw[:Ω]
    end

    ## Sample total uptake capacity per consumer
    if !haskey(kw, :T)
        T = fill(1.0, N)
    else
        T = kw[:T]
    end

    ## Generate uptake matrix from a dirichlet distribution
    for i = 1:N
        dD = Dirichlet(Ω[i] * θ[i, :])
        u[i, :] = rand(dD) * T[i]
    end
    return u

end



function MiCRM_par_therm(N, M; f_u=def_u, f_l=def_l, f_ρ=def_ρ, f_ω=def_ω, f_m=def_m, kwargs...)
    ## We initialize the parameters: uptake matrix, leakage tensor, inflow, and
    ## maitenance biomass

    kw = Dict{Symbol, Any}(kwargs)

    if !haskey(kw, :d)
        dissD = Uniform(0.0, 0.05)
        d = rand(dissD, N, M)
    else
        d = kw[:d]
    end

    ρ = f_ρ(N, M, kw)
    ω = f_ω(N, M, kw)
    m = f_m(N, M, kw)
    l = f_l(N, M, d, kw)
    u = f_u(N, M, kw)

    λ = zeros(N, M)

     ## Calculating total leakage of consumers per resource
    for i in 1:N
        for α in 1:M
            λ[i, α] = sum(l[i, α, :])
        end
    end

    if !haskey(kw, :d)
        kw_nt = (; kwargs...)
        p_nt = (N=N, M=M, u=u, m=m, l=l, ρ=ρ, ω=ω, λ=λ, d=d)
    else
        kw_nt = (; kwargs...)
        p_nt = (N=N, M=M, u=u, m=m, l=l, ρ=ρ, ω=ω, λ=λ)
    end

    out = Base.merge(p_nt, kw_nt)

    return out 

end
