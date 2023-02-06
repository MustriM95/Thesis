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

function CoOp(; p)
    eff_L = zeros(N, M)

    for j in 1:N
        for i in 1:M
            eff_L[j, :] += p.u[j, i]*p.l[j, i, :]
        end
    end

    cfeed(x, y) = dot(x, y) / (norm(x))

    CO = 0.0
    for i in 1:N
        for j in (i+1):N
            CO += cfeed(p.u[i, :], eff_L[j, :])*(2/(N*(N-1)))
        end
    end

    return CO/p.L
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
    ortho = nullspace(θ_uni')

    for j in 1:(N-1)
        for i in 1:M
            if ortho[i, j] < 0.0
                ortho[i, j] = 0.0
            end
        end
    end
    
    null_ind = 1
    min = 0.0
    for i in 1:(N-1)
        dist = dot(ortho[:, i], θ_uni)
        if ==(i, 1)
           min = dist
        else 
            if dist < min
                min = dist
                null_ind = i
            end
        end
    end
    
    null_vec = ortho[:, null_ind]
    
    null_vec = null_vec/norm(null_vec)
    
    ϕ = fill(0.0, M)
    ϕ = θ_uni*sim + null_vec
    ϕ = ϕ/norm(ϕ)

    dD = Dirichlet(100*ϕ[:])

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

p = generate_params(N, M; f_u=F_u, f_m=F_m, f_ρ=F_ρ, f_ω=F_ω, f_l=F_l, L=0.1, θ=θ, Ω = Ω, sim = 100.0)

NiOv(p=p)
CoOp(p=p)

θ_uni = zeros(N)
for i in 1:N
    θ_uni += θ[i, :]/N
end
θ_uni = θ_uni/norm(θ_uni)
ortho =nullspace(θ_uni')

for j in 1:(N-1)
    for i in 1:M
        if ortho[i, j] < 0.0
            ortho[i, j] = 0.0
        end
    end
end

ortho = ortho ./ norm.(eachcol(ortho))'

ortho

null_vec = null_vec/norm(null_vec)

dist_sum = norm(ortho[:, 1] - θ_uni) + norm(ortho[:, 2] - θ_uni)
norm(ortho[:, 1] - θ_uni)/dist_sum
norm(ortho[:, 2] - θ_uni)/dist_sum

ϕ = fill(0.0, M)
ϕ = θ_uni*0.00001 + (norm(ortho[:, 1])/dist_sum)*ortho[:, 1] + (norm(ortho[:, 2])/dist_sum)*ortho[:, 2]
ϕ = θ_uni*0.00001 + ortho[:, 1] 
ϕ = ϕ/norm(ϕ)

dot(ϕ,θ_uni)

ϕ
θ_uni
