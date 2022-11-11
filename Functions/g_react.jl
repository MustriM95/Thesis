"""
Function for g-reactivity calculation given a perturbation transformation (C), Jacobian matrix (J) evaluated
at the steady state solution in question.
"""

function g_react(C, J)

    TransJ = transpose(C)*C*J

    H = Hermitian(TransJ)

    H_ev = eigenvals(H)

    g_r = Eig[maximum((real(H_ev[i]), i) for i in 1:dim)[2]]

    return g_r

end
