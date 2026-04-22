"""
Per-period dissolution probability at state (τ, xg[i]).
Combines cheating-triggered divorce with exogenous termination.
divorce_prob(τ, xgi, sol, param)
"""
function divorce_prob(τ, i, sol, param)
    s = sol.s_star_grid[τ + 1, i]
    VS = sol.VS

    if isapprox(s, 0.0; atol=1e-8) #exog only
        return param.λ
    end

    x = sol.xg[i]
    VM_next = (τ < param.T) ? sol.VM[τ + 2, i] : sol.VM[τ + 1, i]
    W = W_cont(VM_next, VS, param)

    x_cut = find_x_cut(W, sol.VM0_vals, sol.xg, param)
    pc = compute_pc(x_cut, param)

    # Faithful-initiated divorce indicator
    faithful_exit = -param.δ_F + param.β * VS
    ind_fd = (faithful_exit >= W) ? 1.0 : 0.0

    if isapprox(s, 1.0; atol=1e-8)
        # (C,C): endogenous divorce when at least one finds outside match
        endo_divorce = 1.0 - (1.0 - pc)^2
    else
        # Mixed: probability s that each cheats
        prob_both_F = (1.0 - s)^2
        prob_one_C = 2.0 * s * (1.0 - s)
        prob_both_C = s^2

        surv_FF = 1.0
        surv_CF = (1.0 - pc) * (1.0 - ind_fd)
        surv_CC = (1.0 - pc)^2

        endo_divorce = 1.0 - (prob_both_F * surv_FF +
                               prob_one_C * surv_CF +
                               prob_both_C * surv_CC)
    end

    # Total dissolution: combine endogenous and exogenous
    # P(dissolve) = 1 - (1-λ)(1-endo)
    return 1.0 - (1.0 - param.λ) * (1.0 - clamp(endo_divorce, 0.0, 1.0))
end

"""
Compute the steady-state distribution of marriages across states (τ, x).
Forward recursion with flow balance at τ=0.

Returns aggregate statistics: infidelity rate, divorce rate, average marriage duration.
"""
function steady_state_distribution(sol, param)

    xg = sol.xg
    Nx = length(xg)
    dx = xg[2] - xg[1]
    VS = sol.VS

    # mass of matched pairs at state
    n = zeros(param.T + 1, Nx)

    # Inflow
    n_S_unnorm = 1.0
    for i in 1:Nx
        if xg[i] >= sol.x_star
            n[1, i] = param.α * n_S_unnorm * f_pdf(param) * dx
        end
    end

    # forward recursion
    for τ in 0:(param.T - 1)
        for i in 1:Nx
            D = divorce_prob(τ, i, sol, param)
            n[τ + 2, i] = (1.0 - D) * n[τ + 1, i]
        end
    end

    n_married = 2.0 * sum(n)
    total = 1.0 + n_married
    n_S = 1.0 / total
    n_pairs = n ./ total

    #Aggregate statistics
    total_married = 2.0 * sum(n_pairs)

    # Infidelity rate
    cheating = 0.0
    for τ in 0:param.T, i in 1:Nx
        s = sol.s_star_grid[τ + 1, i]
        cheating += 2.0 * s * n_pairs[τ + 1, i]
    end
    infidelity_rate = total_married > 0 ? cheating / total_married : 0.0

    # Divorce rate
    dissolving = 0.0
    for τ in 0:param.T, i in 1:Nx
        D = divorce_prob(τ, i, sol, param)
        dissolving += D * n_pairs[τ + 1, i]
    end
    divorce_rate = total_married > 0 ? dissolving / (total_married / 2.0) : 0.0

    # Average marriage duration
    weighted_dur = 0.0
    total_diss = 0.0
    for τ in 0:param.T, i in 1:Nx
        D = divorce_prob(τ, i, sol, param)
        weighted_dur += τ * D * n_pairs[τ + 1, i]
        total_diss += D * n_pairs[τ + 1, i]
    end
    avg_duration = total_diss > 0 ? weighted_dur / total_diss : Inf

    return (
        n_S = n_S,
        n_pairs = n_pairs,
        infidelity_rate = infidelity_rate,
        divorce_rate = divorce_rate,
        avg_duration = avg_duration,
    )
end
