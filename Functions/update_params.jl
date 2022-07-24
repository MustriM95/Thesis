function update_parameters!(sys::ODESystem, p_new)
    @assert all([:l,:ρ,:N,:M,:m,:ω,:u] .∈ Ref(collect(keys(p_new)))) "missing parameters"
    @assert (p_new[:N] + p_new[:M]) ==  length(states(sys)) "system size is different"

    #overwrite defaults
    for (k,v) in p_new
        if k ∉ [:N, :M]
            for I = CartesianIndices(v)
                if length(I) == 1
                    sys.defaults[ModelingToolkit.toparam(Symbolics.variables(k,I[1]))[1]] = v[I]
                elseif length(I) == 2
                    sys.defaults[ModelingToolkit.toparam(Symbolics.variables(k,I[1],I[2]))[1]] = v[I]
                else
                    error("$k has wrong dimensions")
                end
            end
        end
    end
end
