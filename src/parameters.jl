"""
Three specifications for how match value evolves over tenure.

  - LinearGamma:  u(τ, x) = x + γ·τ                
  - ConcaveGamma: u(τ, x) = x + γ_max·ρ·τ/(1+ρ·τ)  
"""
abstract type GammaSpec end


struct LinearGamma <: GammaSpec
    γ::Float64
end
LinearGamma() = LinearGamma(0.08)

"""
    ConcaveGamma(γ_max, ρ)

Concave, diminishing returns to marriage; approaches γ_max asymptotically.
"""
struct ConcaveGamma <: GammaSpec
    γ_max::Float64   
    ρ::Float64       
end
ConcaveGamma() = ConcaveGamma(1.5, 0.15)



"""
    gamma_val(spec, τ) → Float64

Deterministic component added to x at tenure τ.
"""
gamma_val(g::LinearGamma, τ::Int)  = g.γ * τ
gamma_val(g::ConcaveGamma, τ::Int) = g.γ_max * g.ρ * τ/(1.0 + g.ρ * τ)


"""
    ModelParams

All parameters of the marriage/infidelity model.
  - α:  meeting probability for singles
  - α′: meeting probability for cheating agents
  - λ:  exogenous marriage termination rate
  - c:  cost of cheating on the cheater
  - κ:  cost on the faithful partner when the other cheats 
  - δ_A: amicable divorce cost (both pay)
  - δ_C: cheater's divorce cost
  - δ_F: faithful partner's divorce cost 
"""
@with_kw struct ModelParams
    β::Float64    = 0.95     
    α::Float64    = 0.30    
    α_prime::Float64 = 0.30  

    #Singles flow payoff
    b::Float64    = 0.50

    # Match quality distribution: x ~ Uniform(x_lo, x_hi)
    # For Linear/Concave: x is drawn once and fixed; Nx is the state grid for x.
    # For AR1Gamma: x sets u_0 only; the ongoing u grid is AR1Gamma.n_u.
    x_lo::Float64 = 0.20      
    x_hi::Float64 = 2.00      
    Nx::Int       = 80       

    # Match value dynamics
    gamma_spec::GammaSpec = LinearGamma()
    T::Int = 35

    #Infidelity costs
    c::Float64 = 0.10    
    κ::Float64 = 0.15     

    # Divorce costs
    δ_A::Float64  = 1.00     
    δ_C::Float64  = 1.50      
    δ_F::Float64  = 0.50      

    #Exogenous termination 
    λ::Float64    = 0.02      

    # Solver settings
    tol::Float64  = 1e-6
    max_iter::Int = 500
end
