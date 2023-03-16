
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