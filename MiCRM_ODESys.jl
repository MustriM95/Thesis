"""
This function takes a set of parameters and builds
the MiCRM system of equations
"""
function MiCRM(; name, p)
    @assert all([:l, :ρ, :ω, :m, :M, :N, :u] .∈ Ref(collect(keys(p))))

    M = p[:M]
    N = p[:N]

    l = ModelingToolkit.toparam(Symbolics.variables(:l, 1:N, 1:M, 1:M))
    ρ = ModelingToolkit.toparam(Symbolics.variables(:ρ, 1:M))
    ω = ModelingToolkit.toparam(Symbolics.variables(:ω, 1:M))
    m = ModelingToolkit.toparam(Symbolics.variables(:m, 1:N))
    u = ModelingToolkit.toparam(Symbolics.variables(:u, 1:N, 1:M))

    sts = @variables t C[1:N](t) R[1:M](t)
    D = Differential(t)


    eqns = []

    for i in 1:N
        RHS = -C[i] * m[i]
        for α in 1:M
            RHS += C[i] * R[α] * u[i, α]
            for β in 1:M
                RHS += -C[i] * R[α] * u[i, α] * l[i, α, β]
            end
        end
        push!(eqns, D(C[i]) ~ RHS)
    end

    # Build system of equations for resources

    for β in 1:M
        RHS = ρ[β] - R[β]*ω[β]
        for i in 1:N
            RHS += -u[i, β] * R[β] * C[i]
            for α in 1:M
                RHS += l[i, α, β] * u[i, α] * C[i] * R[α]
            end
        end
        push!(eqns, D(R[β]) ~ RHS)
    end

    vars = vcat(l[:], ρ, ω, m, M, N, u[:])
    vals = vcat(p[:l][:], p[:ρ], p[:ω], p[:m], p[:M], p[:N], p[:u][:])

    sys = ODESystem(eqns, t; name, defaults = Dict(vars[i] => vals[i] for i = eachindex(vars)))

    return(sys)
end
