module MarriageModel

using Parameters
using Statistics
using Interpolations
using QuadGK
using Plots
using Random
using Distributions
using Printf
using Statistics

include("parameters.jl")
include("functions.jl")
include("payoffs.jl")
include("nash.jl")
include("backward_induction.jl")
include("outer_loop.jl")
include("steady_state.jl")

export GammaSpec, LinearGamma, ConcaveGamma, AR1Gamma
export gamma_val
export ModelParams
export x_grid, F_cdf, f_pdf, match_value, W_cont
export find_x_cut, compute_pc, integral_VM0_above, find_x_star
export payoff_FF, payoff_CF, payoff_CC
export nash_equilibrium
export solve_terminal, backward_induction
export compute_VS, solve_model
export divorce_prob, steady_state_distribution



end 
