"""
Function for g-reactivity calculation given a perturbation transformation (C), Jacobian matrix (J) evaluated
at the steady state solution in question.
"""

function g_react(; C = nothing, J)

    if  ==(C,nothing)
        N = dim(J)
        C = 1.0*Matrix(I, N, N)
    end

    TransJ = transpose(C)*C*J

    H = Hermitian(TransJ)

    H_ev = eigvals(H)

    g_r = H_ev[maximum((real(H_ev[i]), i) for i in 1:dim(H))[2]]

    return g_r

end
