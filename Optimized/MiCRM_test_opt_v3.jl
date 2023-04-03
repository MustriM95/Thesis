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
            NV += ((NO - cossim(p.u[i, :], p.u[j, :])^2)*(2/(N*(N-1))))
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
            CV += ((CO - cossim(p.u[i, :], eff_L[j, :])^2)*(2/(N*(N-1))))
        end
    end

    return CV
end

condition(du, t, integrator) = norm(integrator(t, Val{1})) <= 1e-5

affect!(integrator) = terminate!(integrator)

cb = DiscreteCallback(condition, affect!)

function MiC_test(; p, t_span = 10000)

## Calculating niche overlap as average pairwise cosine similarity
    NO = NiOv(p=p)
    NV = NiVar(NO=NO, p=p)
## We calculate the cross feeding/cooperation index
    CO = CoOp(p=p)
    CV = CoVar(CO=CO, p=p)

## Building the ODEs and initializing

    CRi = fill(0.0, (N+M))
    for i in 1:N
        CRi[i] = 0.1
    end

    for α in (N+1):(N+M)
        CRi[α] = 0.1
    end

    prob = ODEProblem(dxx!, CRi, (0.0, t_span), p)

    sol =solve(prob, reltol=1e-9,abstol=1e-9, CVODE_BDF(), callback=cb)

    Jac = MiCRM_jac(p=p, sol=sol)

    gR = g_react(J=Jac)

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

    Ci = fill(0.0, N)
    for i in 1:N
        Ci[i] = 0.1
    end

    prob_LV = ODEProblem(LV_dx!, Ci, (0.0, t_span), LV1)

    sol_LV = solve(prob_LV, reltol=1e-9,abstol=1e-9,
     CVODE_BDF(max_convergence_failures = 3, max_hnil_warns = 1,
      stability_limit_detect=true), callback=cb)

    eq_t = floor(Int, sol.t[end])

    if eq_t == t_span
        println("Steady state not reached, increase t_span")
    end
    
    #while norm(sol(eq_t, Val{1})) > 1e-5 
     #   eq_t += 1
      #  if eq_t == t_span
       #     println("Steady state not reached, increase t_span")
        #    break
        #end
    #end
    

    dtmin_err = sol_LV.retcode
    
    
    if eq_t > sol_LV.t[length(sol_LV)]
        eq_t = floor(Int, sol_LV.t[length(sol_LV)])
    end
    
    println(eq_t)

    C_eq = zeros(N)
    C_LV_eq = zeros(N)

    for i in 1:N
        C_eq[i] = sol(eq_t, idxs = i)
    end

    for i in 1:N
        C_LV_eq[i] = sol_LV(eq_t, idxs = i)
    end

    LV_Jac = Eff_Lv_Jac(p_lv = LV1, sol = sol)

    Eig = eigvals(Jac)

    EigLV = eigvals(LV_Jac)

    gR_LV = g_react(J=LV_Jac)

    domEig = Eig[maximum((real(Eig[i]), i) for i in 1:(N+M))[2]]
    domEigLV = EigLV[maximum((real(EigLV[i]), i) for i in 1:N)[2]]

    SMAPE = 0.0
    for i in 1:N
        for t in 1:Int(eq_t)
            SMAPE += log(abs(sol_LV(t, idxs=i)/sol(t, idxs=i)))/(eq_t*N)
        end
    end


    Eq_SMAPE = 0.0

    for i in 1:N
        Eq_SMAPE += log(abs(C_LV_eq[i]/C_eq[i]))/N
    end

    C_sur = 0.0
    for i in 1:N
        if sol(eq_t, idxs=i) > 1e-6
            C_sur += 1/N
        end
    end


    Sim_res = Dict(:NO => NO, :NV => NV, :CO => CO, :CV => CV, :eq_t => eq_t,
     :domEig => domEig, :domEigLV => domEigLV, :SMAPE => SMAPE, :Eq_SMAPE => Eq_SMAPE,
     :C_sur => C_sur, :trc_max => trc_max, :gR => gR, :gR_LV => gR_LV, :dtmin_err => dtmin_err)
 end
