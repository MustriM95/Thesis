

function MiC_test(μ, σ, L, θ, Ω, N, M)

    NoiseD = Normal(μ, σ)
    NoiseM = rand(NoiseD, N, M)

    t_span = 500.0

    θ = θ + NoiseM

    @named p = MiCRM_par(N = N, M = M, L=L, θ = θ)

    CM = cor(p[:u])

    @named sys1 = MiCRM(p = p)

    CRi = fill(states(sys1)[1] => 0.0, (N+M))
    for i in 1:N
        CRi[i] = states(sys1)[i] => 0.1
    end

    for α in (N+1):(N+M)
        CRi[α] = states(sys1)[α] => 0.1
    end

    prob = ODEProblem(sys1, CRi, (0.0, t_span), [], jac = true)

    sol =solve(prob,reltol=1e-8,abstol=1e-8, saveat=1)

    @named Jac = MiCRM_jac(p=p, symbolic = false, sol=sol)

    @named LV1 = Eff_LV_params(p = p, sol = sol, verbose = true)

    LV1[:ℵ]
    LV1[:r]

    @named LVM = EFF_LV_sys(p_lv = LV1)

    Ci = fill(0.1, N)
    Ci = [states(LVM)[i] => u0[i] for i = eachindex(Ci)]

    prob_LV = ODEProblem(LVM, u0, (0.0, t_span), jac = true)

    sol_LV = solve(prob_LV, reltol=1e-8,abstol=1e-8, saveat=1)

    C_eq = zeros(N)
    C_LV_eq = zeros(N)

    for i in 1:N
        C_eq[i] = sol[i, 300]
    end

    for i in 1:N
        C_LV_eq[i] = sol_LV[i, 300]
    end

    @named LV_Jac = Eff_Lv_Jac(p_lv = LV1, sol = C_eq)

    Eig = eigvals(Jac)

    EigLV = eigvals(LV_Jac)

    domEig = maximum(real(Eig))
    domEigLV = maximum(real(EigLV))

    logL_p = 0.0

    for j in 1:N
        av_sum = 0.0
        for i in 1:length(sol)
            av_sum += abs(sol[j, i] - sol_LV[j, i])/(sol[j, i])
        end
        logL_p = av_sum/N
    end

    logL_null = 0.0

    for j in 1:N
        av_sum = 0.0
        for i in 1:length(sol)
            av_sum += abs(sol[j, i] - sol[j, 0])/(sol[j, i])
        end
        logL_null += av_sum/N
    end


    p_R2 = 1 - (log(logL_p)/log(logL_null))

    Eq_chi = 0.0

    for i in 1:N
        Eq_chi += ((C_eq[i] - C_LV_eq[i])^2)
    end

    Eq_chi
    
