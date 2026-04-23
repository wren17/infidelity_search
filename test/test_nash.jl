using Test
using MarriageModel

p = ModelParams()
xg = x_grid(p)
Nx = length(xg)

VS_test = p.b / (1.0 - p.β)
VM0_vals = [match_value(0, x, p) / (1.0 - p.β) for x in xg]

W_high = VS_test * 3.0
W_low  = VS_test * 0.3

x_mid = (p.x_lo + p.x_hi) / 2.0
τ = 5

@testset "Payoffs" begin

    @testset "payoff_CF" begin
        V1, V2 = payoff_CF(τ, x_mid, W_high, VM0_vals, xg, VS_test, p)

        u = match_value(τ, x_mid, p)
        @test V1 < u + W_high + 1.0   # cheater not getting something for nothing
        @test V2 < payoff_FF(τ, x_mid, W_high, p)  # faithful is hurt by partner cheating

    end

    @testset "payoff_CC" begin
        V = payoff_CC(τ, x_mid, W_high, VM0_vals, xg, VS_test, p)
        @test isfinite(V)
        # CC payoff < FF payoff when W is high
        V_FF = payoff_FF(τ, x_mid, W_high, p)
        @test V < V_FF
    end

    @testset "high W has payoff ordering: V_FF > V_CC" begin
        W_very_high = VS_test * 3.0
        V_FF = payoff_FF(τ, x_mid, W_very_high, p)
        V_CC = payoff_CC(τ, x_mid, W_very_high, VM0_vals, xg, VS_test, p)
        @test V_FF > V_CC
    end

end

@testset "Nash equilibrium" begin

    @testset "FF equilibrium when W is very high" begin
        W_very_high = VS_test * 5.0
        eq_type, s_star, V_eq = nash_equilibrium(τ, x_mid, W_very_high, VM0_vals, xg, VS_test, p)
        @test eq_type in (:FF, :multiple)
        @test isfinite(V_eq)
    end


    @testset "s_star decreasing in x" begin
        x_vals = xg[10:10:end]
        s_vals = [nash_equilibrium(τ, x, W_high * 0.8, VM0_vals, xg, VS_test, p)[2] for x in x_vals]
        @test all(diff(s_vals) .<= 1e-6)
    end


end
