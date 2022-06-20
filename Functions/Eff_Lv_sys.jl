"""
This function takes a dictionary of parameters constructed by
Eff_LV_params and creates a symbolic ODESystem for the LV model

"""

function EFF_LV_sys(; name, p_lv)
    @assert all([:ℵ, :K, :r, :N] .∈ Ref(collect(keys(p_lv))))

    N = p_lv[:N]
    ℵ = ModelingToolkit.toparam(Symbolics.variables(:ℵ, 1:N, 1:N))
    K = ModelingToolkit.toparam(Symbolics.variables(:K, 1:N))
    r = ModelingToolkit.toparam(Symbolics.variables(:r, 1:N))

    sts = @variables t C[1:N](t)

    D = Differential(t)

    eqns = []

    for i in 1:N
        RHS = r[i]*C[i] - (r[i]/K[i])*C[i]*C[i]
        for j in 1:N
            if j != i
                RHS += -r[i]/K[i]*ℵ[i, j]*C[j]*C[i]
            end
        end
        push!(eqns, D(C[i]) ~ RHS)
    end

    vars = vcat(ℵ[:], K, r)
    vals = vcat(p_lv[:ℵ][:], p_lv[:K], p_lv[:r])

    sys = ODESystem(eqns, t; name, defaults = Dict(vars[i] => vals[i] for i = eachindex(vars)))

    return(sys)
end
