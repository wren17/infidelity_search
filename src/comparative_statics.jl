"""
Vary divorce cost asymmetry.
"""
function comparative_statics_δ(param; n_points=5)

    δ_sum = param.δ_C + param.δ_F
    δ_C_vals = range(param.δ_A, δ_sum - 0.01, length=n_points)
    results = []

    for δ_C in δ_C_vals
        δ_F = δ_sum - δ_C
        p = ModelParams(param; δ_C=δ_C, δ_F=δ_F)
        sol = solve_model(p)
        ss = steady_state_distribution(sol, p)

        push!(results, (
            δ_C = δ_C,
            δ_F = δ_F,
            asymmetry = δ_C - δ_F,
            VS = sol.VS,
            x_star = sol.x_star,
            infidelity_rate = ss.infidelity_rate,
            divorce_rate = ss.divorce_rate,
            avg_duration = ss.avg_duration
        ))

        @printf("δ_C=%.2f  δ_F=%.2f  I=%.3f  D=%.3f\n",
                δ_C, δ_F, ss.infidelity_rate, ss.divorce_rate)
    end

    return results
end

"""
Vary the linear capital growth rate γ.
"""
function comparative_statics_γ_linear(param; n_points=10)

    γ_vals = range(0.01, 0.20, length=n_points)
    results = []

    for γ in γ_vals
        p = ModelParams(param; gamma_spec=LinearGamma(γ))
        sol = solve_model(p)
        ss = steady_state_distribution(sol, p)

        push!(results, (
            γ = γ,
            spec = :linear,
            VS = sol.VS,
            x_star = sol.x_star,
            infidelity_rate = ss.infidelity_rate,
            divorce_rate = ss.divorce_rate,
            avg_duration = ss.avg_duration
        ))

        @printf("Linear γ=%.3f  I=%.3f  D=%.3f  dur=%.1f\n",
                γ, ss.infidelity_rate, ss.divorce_rate, ss.avg_duration)
    end

    return results
end

"""
Vary the curvature parameter ρ of ConcaveGamma, holding γ_max fixed.
"""
function comparative_statics_γ_concave(param; n_points=5, γ_max=ConcaveGamma().γ_max)

    ρ_vals = range(0.05, 0.50, length=n_points)
    results = []

    for ρ in ρ_vals
        p = ModelParams(param; gamma_spec=ConcaveGamma(γ_max, ρ))
        sol = solve_model(p)
        ss = steady_state_distribution(sol, p)

        push!(results, (
            ρ = ρ,
            γ_max = γ_max,
            spec = :concave,
            VS = sol.VS,
            x_star = sol.x_star,
            infidelity_rate = ss.infidelity_rate,
            divorce_rate = ss.divorce_rate,
            avg_duration = ss.avg_duration
        ))

        @printf("Concave ρ=%.3f  I=%.3f  D=%.3f  dur=%.1f\n",
                ρ, ss.infidelity_rate, ss.divorce_rate, ss.avg_duration)
    end

    return results
end


"""
Run linear and concave γ specifications side by side.
Returns a Dict with keys :linear, :concave.
"""
function comparative_statics_γ_all(param; n_points=8)
    return Dict(
        :linear => comparative_statics_γ_linear(param; n_points),
        :concave => comparative_statics_γ_concave(param; n_points),
    )
end

"""
Solve the model for several linear γ values
"""
function gamma_heatmap_solutions_linear(param; γ_vals=[0.02, 0.05, 0.08, 0.12])
    sols = []
    params = []
    for γ in γ_vals
        p = ModelParams(param; gamma_spec=LinearGamma(γ))
        sol = solve_model(p)
        push!(sols, sol)
        push!(params, p)
    end
    return sols, params, γ_vals
end

"""
Solve the model for several concave γ ρ values
"""
function gamma_heatmap_solutions_concave(param; ρ_vals=[0.05, 0.10, 0.15, 0.30])
    sols = []
    params = []
    for ρ in ρ_vals
        γ_max = ConcaveGamma().γ_max
        p = ModelParams(param; gamma_spec=ConcaveGamma(γ_max, ρ))
        sol = solve_model(p)
        push!(sols, sol)
        push!(params, p)
    end
    return sols, params, ρ_vals
end

"""
plot_comparative_statics_δ(results)

Plot infidelity and divorce against cost asymmetry.
"""
function plot_comparative_statics_δ(results)

    asym = [r.asymmetry for r in results]
    p1 = plot(asym, [r.infidelity_rate for r in results],
              title="Infidelity vs. δ_C - δ_F", xlabel="δ_C - δ_F",
              ylabel="Infidelity Rate", color=:red, lw=2, marker=:circle, legend=false)
    p2 = plot(asym, [r.divorce_rate for r in results],
              title="Divorce vs. δ_C - δ_F", xlabel="δ_C - δ_F",
              ylabel="Divorce Rate", color=:blue, lw=2, marker=:circle, legend=false)

    return plot(p1, p2, layout=(1, 2), size=(900, 350))
end

"""
plot_comparative_statics_γ(results_dict)

Plot outcomes across the three γ specifications.
"""
function plot_comparative_statics_γ(results_dict)

    res_lin = results_dict[:linear]
    res_con = results_dict[:concave]

    p1 = plot(title="Infidelity Rate by γ Specification",
              xlabel="Parameter value", ylabel="Infidelity Rate",
              legend=:topright, size=(700, 400))
    plot!(p1, [r.γ for r in res_lin], [r.infidelity_rate for r in res_lin],
          label="Linear (γ)", color=:blue, lw=2, marker=:circle)
    plot!(p1, [r.ρ for r in res_con], [r.infidelity_rate for r in res_con],
          label="Concave (ρ)", color=:red, lw=2, marker=:square)


    p2 = plot(title="Divorce Rate by γ Specification",
              xlabel="Parameter value", ylabel="Divorce Rate",
              legend=:topright, size=(700, 400))
    plot!(p2, [r.γ for r in res_lin], [r.divorce_rate for r in res_lin],
          label="Linear (γ)", color=:blue, lw=2, marker=:circle)
    plot!(p2, [r.ρ for r in res_con], [r.divorce_rate for r in res_con],
          label="Concave (ρ)", color=:red, lw=2, marker=:square)

    p3 = plot(title="Average Marriage Duration",
              xlabel="Parameter value", ylabel="Duration",
              legend=:topright, size=(700, 400))
    plot!(p3, [r.γ for r in res_lin], [r.avg_duration for r in res_lin],
          label="Linear (γ)", color=:blue, lw=2, marker=:circle)
    plot!(p3, [r.ρ for r in res_con], [r.avg_duration for r in res_con],
          label="Concave (ρ)", color=:red, lw=2, marker=:square)

    return plot(p1, p2, p3, layout=(1, 3), size=(1500, 400))
end

"""
Heatmap grid comparing s*(τ, x) across different linear γ values.
"""
function plot_gamma_heatmaps_linear(sols, params, γ_vals)
    n = length(sols)
    plts = [heatmap(sols[i].xg, 0:params[i].T, sols[i].s_star_grid,
                    title="Linear γ=$(γ_vals[i])",
                    xlabel="x", ylabel="τ",
                    color=:heat, clims=(0, 1), titlefontsize=10)
            for i in 1:n]
    return plot(plts..., layout=(1, n), size=(300*n, 350))
end

"""
Heatmap grid comparing s*(τ, x) across different concave γ ρ values.
"""
function plot_gamma_heatmaps_concave(sols, params, ρ_vals)
    n = length(sols)
    plts = [heatmap(sols[i].xg, 0:params[i].T, sols[i].s_star_grid,
                    title="Concave ρ=$(ρ_vals[i])",
                    xlabel="x", ylabel="τ",
                    color=:heat, clims=(0, 1), titlefontsize=10)
            for i in 1:n]
    return plot(plts..., layout=(1, n), size=(300*n, 350))
end

"""
Solve the model for several concave γ_max values (holding ρ fixed).
"""
function gamma_heatmap_solutions_concave_γmax(param; γmax_vals=[0.5, 2.0, 5.0, 15.0])
    ρ = ConcaveGamma().ρ
    sols = []
    params = []
    for γ_max in γmax_vals
        p = ModelParams(param; gamma_spec=ConcaveGamma(γ_max, ρ))
        sol = solve_model(p)
        push!(sols, sol)
        push!(params, p)
    end
    return sols, params, γmax_vals
end

"""
Heatmap grid comparing s*(τ, x) across different concave γ_max values.
"""
function plot_gamma_heatmaps_concave_γmax(sols, params, γmax_vals)
    n = length(sols)
    plts = [heatmap(sols[i].xg, 0:params[i].T, sols[i].s_star_grid,
                    title="Concave γ_max=$(γmax_vals[i])",
                    xlabel="x", ylabel="τ",
                    color=:heat, clims=(0, 1), titlefontsize=10)
            for i in 1:n]
    return plot(plts..., layout=(1, n), size=(300*n, 350))
end




