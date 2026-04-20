"""
Payoffs of both players faithful
"""
function payoff_FF(τ, x, W, param)
    return match_value(τ, x, param) + W
end

"""
(C, F)
V1_CF (cheater)
V2_CF (faithful)
"""
function payoff_CF(τ, x, W, VM0_vals, xg, VS, param)

    u = match_value(τ, x, param)

    # cheater reservation quality and probability of finding outside match
    x_cut = find_x_cut(W, VM0_vals, xg, param)
    pc = compute_pc(x_cut, param)


    integral_above = integral_VM0_above(x_cut, VM0_vals, xg, param)
    prob_above = 1.0 - F_cdf(x_cut, param)
    cheater_leave_gain = param.α_prime * (param.β * integral_above - param.δ_C * (1.0 - F_cdf(x_cut, param)))


    # Faithful divorces when
    faithful_exit = -param.δ_F + param.β * VS
    ind_fd = (faithful_exit >= W) ? 1.0 : 0.0

    V1_CF = u - param.c + cheater_leave_gain +
            (1.0 - pc) * (ind_fd * (-param.δ_C + param.β * VS) +
                          (1.0 - ind_fd) * W)

    # faithful payoff
    opt_stay = pc * faithful_exit + (1.0 - pc) * W
    opt_leave = faithful_exit

    V2_CF = u - param.κ + max(opt_leave, opt_stay)

    return V1_CF, V2_CF
end

"""
Payoffs of both players cheating
"""
function payoff_CC(τ, x, W, VM0_vals, xg, VS, param)

    u = match_value(τ, x, param)

    x_cut = find_x_cut(W, VM0_vals, xg, param)
    pc    = compute_pc(x_cut, param)

    integral_above = integral_VM0_above(x_cut, VM0_vals, xg, param)
    prob_above = 1.0 - F_cdf(x_cut, param)
    E_VM0_cond = prob_above > 1e-12 ? integral_above / prob_above : 0.0

    V_CC = u - param.c +
           pc * (-param.δ_C + param.β * E_VM0_cond) +
           (1.0 - pc) * pc * (-param.δ_C + param.β * VS) +
           (1.0 - pc)^2 * W

    return V_CC
end
