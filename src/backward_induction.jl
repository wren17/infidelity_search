"""
After period T, the marriage is stationary: match quality stops improving. Value of marriage stays the smame for all τ ≥ T.

  1. Start with a guess V (closed-form value under forever-faithful). 
  This should be the case if T is large enough to make cheating unattractive.
  2. Compute continuation value W for the foever faithful after T outcome
  3. Solve the Nash equilibrium given W and get the equilibrium payoff V_new.
  4. If V_new is approx. V, done. Otherwise update V and repeat.

If the equilibrium is faithful (FF) for all x, the initial guess is already the fixed point and the loop exits on the first iteration.
"""
function solve_terminal(VM0_vals, xg, VS, param)

    Nx = length(xg)
    VM_T = zeros(Nx)
    eq_types = Vector{Symbol}(undef, Nx)
    s_stars = zeros(Nx)

    for i in 1:Nx
        x = xg[i]

        #assuming that T is big enough to make cheating unattractive
        eq_type = :FF
        s_star = 0.0
        u_base = match_value(param.T, x, param)
        V = (u_base + param.β * param.λ * VS)/(1.0 - param.β * (1.0 - param.λ)) #closeform: geometric sum of forever faithful

        for _ in 1:param.max_iter
            W = W_cont(V, VS, param)
            eq_type, s_star, V_new = nash_equilibrium(param.T, x, W, VM0_vals, xg, VS, param)

            if abs(V_new - V) < param.tol
                V = V_new
                break
            end
            V = V_new #add damping maybe?
        end

        VM_T[i] = V
        eq_types[i] = eq_type
        s_stars[i] = s_star
    end

    return VM_T, eq_types, s_stars
end

"""
Compute V_M(τ, x) for all τ < T and x
via backward induction from the terminal condition.
"""
function backward_induction(VM0_vals, xg, VS, param)

    Nx = length(xg)
    NT = param.T + 1

    VM = zeros(NT, Nx)
    eq_type_grid = Matrix{Symbol}(undef, NT, Nx)
    s_star_grid = zeros(NT, Nx)

    # terminal solution
    VM_T, eq_T, s_T = solve_terminal(VM0_vals, xg, VS, param)
    VM[NT, :] = VM_T
    eq_type_grid[NT, :] = eq_T
    s_star_grid[NT, :] = s_T

    # Backward induction
    for τ in (param.T - 1):-1:0
        for i in 1:Nx
            x = xg[i]
            W = W_cont(VM[τ + 2, i], VS, param)
            eq_type, s_star, V_new = nash_equilibrium(τ, x, W, VM0_vals, xg, VS, param)
            VM[τ + 1, i] = V_new
            eq_type_grid[τ + 1, i] = eq_type
            s_star_grid[τ + 1, i] = s_star
        end
    end

    return VM, eq_type_grid, s_star_grid
end
