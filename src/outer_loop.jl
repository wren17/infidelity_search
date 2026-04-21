
"""
Update V_S using the Bellman equation
"""
function compute_VS(VM0_vals, xg, VS_old, param)

    density = f_pdf(param)
    itp = linear_interpolation(xg, VM0_vals, extrapolation_bc=Flat())
    x_star = find_x_star(VM0_vals, VS_old, xg, param)

    integral, _ = quadgk(x -> itp(x) * density, x_star, param.x_hi; rtol=1e-6)

    F_xstar = F_cdf(x_star, param)
    numer = param.b + param.β * param.α * integral
    denom = 1.0 - param.β + param.β * param.α * (1.0 - F_xstar)

    return numer/denom
end

"""
Outer fixed-point iteration over V_S,
with inner backward induction at each step.
"""
function solve_model(param)

    xg = x_grid(param)
    Nx = length(xg)

    # Initialization
    VS = param.b/(1.0 - param.β)
    VM0_vals = [match_value(0, x, param) / (1.0 - param.β) for x in xg]

    VM = zeros(param.T + 1, Nx)
    eq_type_grid = Matrix{Symbol}(undef, param.T + 1, Nx)
    s_star_grid = zeros(param.T + 1, Nx)

    converged = false

    for iter in 1:param.max_iter

        # backward induction
        VM, eq_type_grid, s_star_grid = backward_induction(VM0_vals, xg, VS, param)
        new_VM0_vals = VM[1, :]

        # Update V_S
        new_VS = compute_VS(new_VM0_vals, xg, VS, param)

        # convergence check
        err_V = maximum(abs.(new_VM0_vals .- VM0_vals))
        err_VS = abs(new_VS - VS)

        if err_V < param.tol && err_VS < param.tol
            converged = true
            VS = new_VS
            VM0_vals = new_VM0_vals
            if iter > 1
                println("Converged in $iter outer iterations.")
            end
            break
        end

        # Dampened update for stability
        VM0_vals = 0.5 .* VM0_vals .+ 0.5 .* new_VM0_vals
        VS = 0.5 * VS + 0.5 * new_VS
    end

    if !converged
        @warn "Outer fixed point did not converge after $(param.max_iter) iterations."
    end

    # Reservation quality
    x_star = find_x_star(VM0_vals, VS, xg, param)

    return (
        VS = VS,
        VM = VM,
        VM0_vals = VM0_vals,
        eq_type_grid = eq_type_grid,
        s_star_grid = s_star_grid,
        x_star = x_star,
        xg = xg,
        converged = converged
    )
end
