"""
This function calculates the Effective Lotka-Volterra
System of ODEs

"""

function EFF_LV_sys(; name, p_lv)
    @assert all([:ℵ, :r, :N] .∈ Ref(collect(keys(p_lv))))

    N = p_lv[:N]
    ℵ = ModelingToolkit.toparam(Symbolics.variables(:ℵ, 1:N, 1:N))
    r = ModelingToolkit.toparam(Symbolics.variables(:r, 1:N))

    sts = @variables t C[1:N](t)

    D = Differential(t)

    eqns = []

    for i in 1:N
        RHS = r[i]*C[i]
        for j in 1:N
            RHS += C[i]*ℵ[i, j]*C[j]
        end
        push!(eqns, D(C[i]) ~ RHS)
    end

    vars = vcat(ℵ[:], r)
    vals = vcat(p_lv[:ℵ][:], p_lv[:r])

    sys = ODESystem(eqns, t; name, defaults = Dict(vars[i] => vals[i] for i = eachindex(vars)))

    return(sys)
end
