# A Dynamic Search Model for Marriage and Infidelity

*Wuyang Ren*

The project extends Burdett, Imai, and Wright (2002) into a discrete-time bilateral strategic model where married agents choose each period whether to remain faithful or search for outside matches. The model also builds on Jovanovic (1979) by letting the flow utility of marriage grow over time, similar to "job-specific" human capital accumulation.

As a course project, the goal is to practice implementations of search models, solving dynamic programming problems through backward induction, and simple Nash equilibrium analysis using computational tools. Some of the tools used here include fixed-point iteration, function approximation, and numerical integration. This framework be relevant for various topics of future research, especially in the field of labour economics.

The model is preliminary; feedbacks and suggestions are welcome.

## The Model

A continuum of agents are either single or married. Single agents receive flow payoff $b$ and meet potential partners at rate $\alpha$, drawing match quality $x \sim F$. Married couples at state $(\tau, x)$ receive flow utility $u(\tau, x) = x + \gamma(\tau)$ and simultaneously choose to be faithful (F) or cheat (C). The flow utility of marriage $u$ grows in $\tau$ for $\tau \leq T$. Cheating costs $c$ to the cheater and imposes $\kappa$ on the faithful partner. Outside matches that would want to match with the cheater arrive at rate $\alpha' \leq \alpha$. Divorce costs are asymmetric, where the cheater bears a higher cost $\delta_C \geq \delta_F$. Marriages also face exogenous termination at rate $\lambda$.

## Solution Algorithm

1.  Discretize match values $x$ on a grid.
2.  Starts with an initial guess for the value of beging single $V_S$ and the value of a new relationship with match value $V_m(0, x)$ .
3.  Solve the cheating game for married couples in each state $(\tau, x)$ using backward induction from period $T$.
4.  Use the solution from backward induction to update $V_S$ and $V_M(0, x)$

## Quick Start

[`InfidelitySearch.qmd`](InfidelitySearch.qmd) provides a more detailed summary of the model and goes through each step of the solution algorithm. If you just want to see the outputs without running anything, an HTML version ([`InfidelitySearch.html`](InfidelitySearch.html)) with the results is also included.

Or, if you'd just like to see some of the model solution, it is also easy to get the results by running the individual components:

```{julia}
using Pkg
Pkg.activate(".")

using MarriageModel

#Pkg.test()

# Baseline (linear γ)
p = ModelParams() #or change the default parameters e.g. ModelParams(Nx=80, T=45)
sol = solve_model(p)
ss = steady_state_distribution(sol, p)

println("\n-- Baseline: Linear γ --")
@printf("V_S = %.4f\n", sol.VS)
@printf("x* = %.4f\n", sol.x_star)
@printf("V_M(0, x_lo) = %.4f\n", sol.VM[1, 1])
@printf("V_M(0, x_hi) = %.4f\n", sol.VM[1, end])

println("\n-- Linear γ Steady State Distributions--")
@printf("Infidelity = %.3f  Divorce = %.3f  Dur = %.3f\n",
         ss.infidelity_rate, ss.divorce_rate, ss.avg_duration)

# Concave γ
p_c = ModelParams(p; gamma_spec=ConcaveGamma(1.5, 0.15))
sol_c = solve_model(p_c)
ss_c = steady_state_distribution(sol_c, p_c)

println("\n-- Concave γ --")
@printf("V_S = %.4f\n", sol_c.VS)
@printf("x* = %.4f\n", sol_c.x_star)
@printf("V_M(0, x_lo) = %.4f\n", sol_c.VM[1, 1])
@printf("V_M(0, x_hi) = %.4f\n", sol_c.VM[1, end])

println("\n-- Concave γ Steady State Distributions--")
@printf("Infidelity = %.3f  Divorce = %.3f  Dur = %.3f\n",
         ss_c.infidelity_rate, ss_c.divorce_rate, ss_c.avg_duration)

```

## Disclosure of AI Use

Claude AI (Sonnet 4.6/Opus 4.6) is used to help generate and fix codes, and revise texts for clarity. It is used more heavily where the tasks are more repetitive: changing the parameters in comparative statics, making multiple graphs etc. All ideas are my own.

## References

-   Burdett, K., Imai, R., & Wright, R. D. (2002). Unstable Relationships. *SSRN Electronic Journal*. https://doi.org/10.2139/ssrn.337201

-   Jovanovic, B. (1979). Job Matching and the Theory of Turnover. *Journal of Political Economy*, *87*(5, Part 1), 972–990. https://doi.org/10.1086/260808