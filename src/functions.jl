"""
Flow utility of a marriage with initial quality x at τ:
"""
function match_value(τ, x, param)
    return x + gamma_val(param.gamma_spec, τ)
end

"""
Discretise match quality grid on [x_lo, x_hi].
"""
function x_grid(param)
    return collect(range(param.x_lo, param.x_hi, length=param.Nx))
end

"""
CDF of uniform(x_lo, x_hi)
"""
function F_cdf(x, param)
    return clamp((x - param.x_lo) / (param.x_hi - param.x_lo), 0.0, 1.0)
end

"""
PDF of uniform(x_lo, x_hi)
"""
f_pdf(param) = 1.0 / (param.x_hi - param.x_lo)

"""
Discounted continuation value of marriage
"""
function W_cont(V_M_next, VS, param)
    return param.β * ((1.0 - param.λ) * V_M_next + param.λ * VS)
end

"""
Find the cheater's reservation quality x_C(τ,x)

Returns x_hi if no solution (cheater never leaves), x_lo if always leaves.
"""
function find_x_cut(W, VM0_vals, xg, param) #VM0_vals: valuation of marriage at τ=0 for each x in xg
    target = (W + param.δ_C) / param.β # cheater's targeted value to justify leaving the relatinship

    # Check boundary
    if target <= VM0_vals[1]
        return param.x_lo    # even lowest x' is worth leaving for
    elseif target >= VM0_vals[end]
        return param.x_hi    # no x' is worth leaving for
    end

    # Linear search + interpolation
    for i in 1:(length(xg)-1)
        if VM0_vals[i] <= target && VM0_vals[i+1] >= target
            frac = (target - VM0_vals[i]) / (VM0_vals[i+1] - VM0_vals[i])
            return xg[i] + frac * (xg[i+1] - xg[i])
        end
    end

    return param.x_hi
end

"""
Probability that a cheating spouse finds a match worth leaving for
"""
function compute_pc(x_cut, param)
    return param.α_prime * (1.0 - F_cdf(x_cut, param))
end

"""
Find the reservation match quality x* where V_M(0, x*) = VS.
"""
function find_x_star(VM0_vals, VS, xg, param)
    target = VS

    if target <= VM0_vals[1]
        return param.x_lo
    elseif target >= VM0_vals[end]
        return param.x_hi
    end

    for i in 1:(length(xg) - 1)
        if VM0_vals[i] <= target && VM0_vals[i + 1] >= target
            frac = (target - VM0_vals[i]) / (VM0_vals[i + 1] - VM0_vals[i])
            return xg[i] + frac * (xg[i + 1] - xg[i])
        end
    end

    return param.x_hi
end

"""
The integral part of the cheater's expected gain from an outside draw
"""
function integral_VM0_above(x_cut, VM0_vals, xg, param)
    a = clamp(x_cut, param.x_lo, param.x_hi)
    itp = linear_interpolation(xg, VM0_vals, extrapolation_bc=Flat())
    density = f_pdf(param)
    I, _ = quadgk(x -> itp(x) * density, a, param.x_hi; rtol=1e-8)
    return I
end
