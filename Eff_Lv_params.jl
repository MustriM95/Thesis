
"""
This function takes MiCRM parameters and
equilibrium solutions and calculates the parameters
for the effective lotka-volterra system. It returns
the parameters in a dictionary of keywords.

"""

function Eff_LV_params(; name, p, sol)
    @assert all([:l, :ρ, :ω, :m, :M, :N, :u] .∈ Ref(collect(keys(p)))) "missing parameters"

    M = p[:M]
    N = p[:N]
    l = p[:l]
    ρ = p[:ρ]
    ω = p[:ω]
    m = p[:m]
    u = p[:u]


    Ceq = sol[1:N, length(sol)]
    Req = sol[(N+1):(N+M), length(sol)]
    A = zeros(M, M)
    ∂R = zeros(M, N)
    ℵ₊ = zeros(N, N)
    ℵ₋ = zeros(N)
    ℵ = zeros(N, N)
    λ = zeros(N, M)
    O₊ = zeros(N) # Dummy variable for numerator of first term of sum
    O = zeros(N) # Dummy variable for first term of sum
    M_lv = zeros(N) # Dummy variable for maitenance term
    P = zeros(N) # Dummy variable for interaction term
    K = zeros(N) # Effective Carrying capacities
    r = zeros(N) # Effective growth rates

## Calculating the A matrix from parameters at equilibrium
## solutions
    for α in 1:M
        for β in 1:M
            A[α, β] = -ω[α]
            for i in 1:N
                A[α, β] += l[i, α, β] * u[i, β] * Ceq[i] -
                u[i, β] * Ceq[i] * δ(α, β)
            end
        end
    end

## Calculating partial derivatives of eq resources with respect
## to consumers
    for α in 1:M
        for j in 1:N
            for β in 1:M
                for γ in 1:M
                    ∂R[α, j] += inv(A)[α, β] * u[j, β] * Req[γ]*(δ(β, γ)
                    - l[j, β, γ])
                end
            end
        end
    end

## Calculating total leakage of consumers per resource
    for i in 1:N
        for α in 1:M
            λ[i, α] = sum(l[i, α, :])
        end
    end

## Calculating numerator terms of interaction matrix
    for i in 1:N
        for j in 1:N
            for α in 1:M
                ℵ₊[i, j] += u[i, α] * (1 - λ[i, α]) * ∂R[α, j]
            end
        end
    end

## Calculating denominator terms of interaction matrix
    for i in 1:N
        for β in 1:M
            ℵ₋[i] += u[i, β] * (1 - λ[i, β]) * ∂R[β, i]
        end
    end

## Calculating interaction matrix
    for i in 1:N
        for j in 1:N
            ℵ[i, j] = ℵ₊[i, j]/ℵ₋[i]
        end
    end

## Calculating effective carrying capacities
    for i in 1:N
        for α in 1:M
            O₊[i] += u[i, α]*(1 - λ[i, α])*Req[α]
        end
    end

    for i in 1:N
        O[i] = O₊[i]/ℵ₋[i]
    end

    for i in 1:N
        M_lv[i] = m[i]/ℵ₋[i]
    end

    for i in 1:N
        for j in 1:N
            P[i] += ℵ[i, j]*Ceq[j]
        end
    end

    for i in 1:N
        K[i] = -O[i] - P[i] + M_lv[i]
    end
## Calculating effective intrinsic growth rates
    for i in 1:N
        r[i] = -K[i]*ℵ₋[i]
    end

    return Dict(:ℵ => ℵ, :K => K, :r => r, :N => N)
end
