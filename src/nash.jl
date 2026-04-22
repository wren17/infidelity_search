"""
Solving the symmetric game at each marriage state (\\tau, x)

Returns:
  - eq_type: :FF, :CC, :mixed, or :multiple
  - s_star:  equilibrium cheating probability (0 for FF, 1 for CC, interior for mixed)
  - V_eq:    equilibrium value for each agent
nash_equilibrium(τ, x, W, VM0_vals, xg, VS, param)  
"""
function nash_equilibrium(τ, x, W, VM0_vals, xg, VS, param)

    V_FF = payoff_FF(τ, x, W, param)

    V1_CF, V2_CF = payoff_CF(τ, x, W, VM0_vals, xg, VS, param)

    V_CC = payoff_CC(τ, x, W, VM0_vals, xg, VS, param)

    V_FC = V2_CF

    #Pure strategy NE conditions:

    # (F,F) is NE when deviating to C is not profitable: V1_CF ≤ V_FF
    ff_is_ne = (V1_CF <= V_FF)

    # (C,C) is NE when deviating to F is not profitable: V_FC ≤ V_CC
    cc_is_ne = (V_FC <= V_CC)

    if ff_is_ne && !cc_is_ne
        return (:FF, 0.0, V_FF)
    elseif cc_is_ne && !ff_is_ne
        return (:CC, 1.0, V_CC)
    elseif ff_is_ne && cc_is_ne
        # Both are NE — select payoff-dominant (higher V)
        if V_FF >= V_CC
            return (:multiple, 0.0, V_FF)
        else
            return (:multiple, 1.0, V_CC)
        end
    end

# Mixed strategy NE: solve for s* in [0,1] such that cheater is indifferent between C and F

    num   = V1_CF - V_FF
    denom = (V1_CF - V_FF) - (V_CC - V_FC)

    if isapprox(denom, 0.0; atol=1e-6)
        # Degenerate default to FF
        return (:FF, 0.0, V_FF)
    end

    s_star = num / denom
    s_star = clamp(s_star, 0.0, 1.0)

    # Equilibrium value at the mixed NE
    V_eq = (1.0 - s_star) * V_FF + s_star * V_FC

    return (:mixed, s_star, V_eq)
end
