"""
Function to generate jacobian for MiCRM system of equations

"""

function MiCRM_jac(; name, p, symbolic = true, sol = nothing)
    @assert all([:l, :ρ, :ω, :m, :M, :N, :u] .∈ Ref(collect(keys(p))))

    M = p[:M]
    N = p[:N]

    l = ModelingToolkit.toparam(Symbolics.variables(:l, 1:N, 1:M, 1:M))
    ρ = ModelingToolkit.toparam(Symbolics.variables(:ρ, 1:M))
    ω = ModelingToolkit.toparam(Symbolics.variables(:ω, 1:M))
    m = ModelingToolkit.toparam(Symbolics.variables(:m, 1:N))
    u = ModelingToolkit.toparam(Symbolics.variables(:u, 1:N, 1:M))

    sts = @variables t C[1:N](t) R[1:M](t)

    if sol != nothing
        C = sol[1:N, length(sol)]
        R = sol[(N+1):(N+M), length(sol)]
    end


    Jac = zeros(Num, (N+M), (N+M))
    λ = zeros(Num, N, M)

    for i in 1:N
        for α in 1:M
            λ[i, α] = sum(l[i, α, :])
        end
    end


    for i in 1:N
        Jac[i, i] = -m[i]
        for α in 1:M
            Jac[i, i] += (1 - λ[i, α])*u[i, α]*R[α]
        end
    end
    for i in 1:N
        for j in (i+1):N
            Jac[i, j] = 0
        end
    end

    for i in 1:N
        for j in (i+1):N
            Jac[j, i] = 0
        end
    end

    for i in 1:N
        for α in (N+1):(N+M)
            Jac[i, α] = C[i]*(1-λ[i, α-N])*u[i, α-N]
        end
    end

    for α in (N+1):(N+M)
        Jac[α, α] = -ω[α-N]
        for i in 1:N
            Jac[α, α] += C[i]*u[i, α-N]*(l[i, α-N, α-N] - 1)
        end
    end

    for α in (N+1):(N+M)
        for i in 1:N
            Jac[α, i] = -u[i, α-N]*R[α-N]
            for β in 1:M
                Jac[α, i] += l[i, β, α-N]*u[i, β]*R[β]
            end
        end
    end

    for α in (N+1):(N+M)
        for β in (N+1):(α-1)
            for i in 1:N
                Jac[α, β] += l[i, β-N, α-N]*u[i, β-N]*C[i]
            end
        end
    end

    for α in (N+1):(N+M)
        for β in (α+1):(N+M)
            for i in 1:N
                Jac[α, β] += l[i, β-N, α-N]*u[i, β-N]*C[i]
            end
        end
    end

    if symbolic == false
        vals= vcat(p[:l][:], p[:ρ], p[:ω], p[:m], p[:u][:])
        vars = vcat(l[:], ρ, ω, m, u[:])

        pmap = [param =>p for (param,p) in zip(vars, vals)]
        Jac = substitute.( Jac, (pmap,))
        Jac = Symbolics.value.(Jac)
    end

    return(Jac)

end
