function MiC_test(; μ, σ, L, N, M, θ = nothing, Ω = nothing, t_span = nothing)

    NoiseD = Normal(μ, σ)
    NoiseM = rand(NoiseD, N, M)
    dU = Uniform(0, 1)

    if t_span == nothing
        t_span = 1000.0
    end

    ## Sample concentration parameters for each consumer from uniform distribution
    if θ == nothing
        θ = zeros(N, M)

        for i in 1:N
            θ[i, :] = rand(dU, M)
        end
    end

    if Ω == nothing
        ## Sample specialisation parameter for each consumer
        Ω = fill(1.0, N)
    end

    θ = θ + broadcast(abs, NoiseM)

    p = generate_params(N, M; f_u=F_u, f_m=F_m, f_ρ=F_ρ, f_ω=F_ω, L=L, θ=θ, Ω = Ω)

## Calculating the average euclidean distances between niches and the niche
# incumbencies
    niche_dist_avg = 0.0

    for i in 1:N
        for j in (i+1):N
            niche_dist_avg += norm(p.u[i, :] - p.u[j, :])*(2/(N*(N-1)))
        end
    end

    NI = log(sqrt(2)/niche_dist_avg)

## Calculating niche overlap as average pairwise cosine similarity
    cossim(x, y) = dot(x, y) / (norm(x)*norm(y))

    NO = 0.0
    for i in 1:N
        for j in (i+1):N
            NO += cossim(p.u[i, :], p.u[j, :])*(2/(N*(N-1)))
        end
    end
## We calculate the cross feeding/cooperation index
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

## Building the ODEs and initializing

    CRi = fill(0.0, (N+M))
    for i in 1:N
        CRi[i] = 0.1
    end

    for α in (N+1):(N+M)
        CRi[α] = 0.1
    end

    prob = ODEProblem(dxx!, CRi, (0.0, t_span), p)

    cb = AutoAbstol()

    sol =solve(prob, reltol=1e-9,abstol=1e-9, saveat=1, TRBDF2(), callback=cb, kwargshandle=KeywordArgSilent)

    Jac = MiCRM_jac(p=p, sol=sol)

    tc = zeros(N)
    tR = zeros(M)

    [tc[i]  = 1/abs(Jac[i, i]) for i in 1:N]
    [tR[α] = 1/abs(Jac[α+N, α+M]) for α in 1:M]

    trc = zeros(N, M)
    for i in 1:N
        for α in 1:M
            trc[i, α] = tR[α]/(tc[i] + tR[α])
        end
    end

    trc_max = maximum(filter(!isnan,trc))


    LV1 = Eff_LV_params(p=p, sol=sol)

    ℵ_m = mean(LV1.ℵ)
    r_m = mean(LV1.r)

    Ci = fill(0.0, N)
    for i in 1:N
        Ci[i] = 0.1
    end

    prob_LV = ODEProblem(LV_dx!, Ci, (0.0, t_span), LV1)

    sol_LV = solve(prob_LV, reltol=1e-9,abstol=1e-9, saveat=1, TRBDF2(), callback=cb, kwargshandle=KeywordArgSilent)

    eq_t = 0
    while norm(sol(eq_t, Val{1})) > 1e-5
        eq_t += 1
        if eq_t == length(sol)
            println("Steady state not reached, increase t_span")
            break
        end
    end

    max_iters = false

    if eq_t > length(sol_LV)
        eq_t = length(sol_LV)
        max_iters = true
    end

    C_eq = zeros(N)
    C_LV_eq = zeros(N)

    for i in 1:N
        C_eq[i] = sol[i, eq_t]
    end

    for i in 1:N
        C_LV_eq[i] = sol_LV[i, eq_t]
    end

    LV_Jac = Eff_Lv_Jac(p_lv = LV1, sol = sol)

    Eig = eigvals(Jac)

    EigLV = eigvals(LV_Jac)

    domEig = maximum(real(Eig))
    domEigLV = maximum(real(EigLV))

    MSE = 0.0
    for i in 1:N
        for t in 1:eq_t
            MSE += (sol[i, t] - sol_LV[i, t])^2/(eq_t*N*0.01)
        end
    end


    Eq_MSE = 0.0

    for i in 1:N
        Eq_MSE += ((C_eq[i] - C_LV_eq[i])^2)/(0.01*N)
    end

    C_sur = 0.0
    for i in 1:N
        if sol[i, eq_t] > 1e-6
            C_sur += 1/N
        end
    end


    Sim_res = Dict(:NO => NO, :CO => CO, :ℵ_m => ℵ_m, :r_m => r_m, :eq_t => eq_t,
     :domEig => domEig, :domEigLV => domEigLV, :MSE => MSE, :Eq_MSE => Eq_MSE,
     :C_sur => C_sur, :trc_max => trc_max)
 end
