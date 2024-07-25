"""

Function to calculate the jacobian of an effective Lotka-Volterra
System. Takes in parameter dictionary from Eff_LV_params and an
optional solution for equilibrium values to find lenarization near
steady state.

"""

function Eff_Lv_Jac(; p_lv, sol)
## Checks that all necessary parameters have been supplied

## Loads in size of system N, and converts the interaction matrix
# and intrinsic growth rates to symbolic variables for MTK.
    N = p_lv.N
    ℵ = p_lv.ℵ
    r = p_lv.r

## We define our state variables
    C = sol[1:N, length(sol)]

    LV_Jac = zeros(N, N)

    LV_Jac = [ℵ[i, j]*C[i] for i in 1:N, j in 1:N]
# reset diagonal
    LV_Jac[diagind(LV_Jac)] .= [r[i] + ℵ[i, i]*C[i] + sum(ℵ[i, j]*C[j] for j in 1:N) for i in 1:N]

    return(LV_Jac)

end
