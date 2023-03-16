#=
Testing niche overlap and cross feeding distributions
=#
cd("C:\\Users\\micho\\github\\Thesis")

using Distributions, LinearAlgebra
include("micrm_params.jl")

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

function NiVar(; NO, p)
    cossim(x, y) = dot(x, y) / (norm(x)*norm(y))

    NV = 0.0
    for i in 1:N
        for j in (i+1):N
            NV += (NO - cossim(p.u[i, :], p.u[j, :])*(2/(N*(N-1))))^2
        end
    end

    return NV
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

function CoVar(; CO, p)
    eff_L = zeros(N, M)

    for j in 1:N
        for i in 1:M
            eff_L[j, :] += p.u[j, i]*p.l[j, i, :]
        end
    end

    cossim(x, y) = dot(x, y) / (norm(x)*norm(y))

    CV = 0.0
    for i in 1:N
        for j in (i+1):N
            CV += (CO - cossim(p.u[i, :], eff_L[j, :])*(2/(N*(N-1))))^2
        end
    end

    return CV
end

############################################################

function F_l(N, M, kw)
    l = zeros(N, M, M)

    if !haskey(kw, :sim)
        sim = 1.0
    else
        sim = kw[:sim]
    end

    if !haskey(kw, :θ)
        dU = Uniform(0, 1)
        ϕ = rand(dU, M)
    else
        θ = kw[:θ]
    end

    θ_uni = zeros(N)
        for i in 1:N
            θ_uni += θ[i, :]/N
        end
    θ_uni
    min_index = 1
        
    for i in 2:N
        temp = θ_uni[1]
        if θ_uni[i] < temp
            min_index = i
        end
    end
        
    ortho = zeros(N)
    ortho[min_index] += θ_uni[min_index]

    u_null = ones(M)

    ϕ = fill(0.0, M)
    #ϕ = θ_uni*sim + (1/sim)*(u_null - θ_uni/norm(θ_uni))
    ϕ = θ_uni*sim + (1/sim)*ortho
    ϕ = ϕ/norm(ϕ)

    dD = Dirichlet(10*ϕ[:])

    for i = 1:N
        for α = 1:M
            l[i, α, :] = rand(dD) * kw[:L]
        end
    end
    return l
end

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
        Ω = kw[:Ω]
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

N=6
M=6

noise = Normal(1.01, 0.1)
NoiseM = rand(noise, N, M)
dU = Uniform(0, 1)

θ = 1.0*I(N)

θ_het = fill(1.0, N, M)

θ = 10.0*θ + broadcast(abs, NoiseM)

Ω = fill(100.0, N)

p = generate_params(N, M; f_u=F_u, f_m=F_m, f_ρ=F_ρ, f_ω=F_ω, f_l=F_l, L=0.5, θ=θ, Ω = Ω, sim = 0.2)

NO = NiOv(p=p)
CO = CoOp(p=p)

NV = NiVar(NO=NO, p=p)
CV = CoVar(CO=CO, p=p)
