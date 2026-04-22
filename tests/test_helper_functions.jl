using Test
using MarriageModel


@testset "x_grid" begin
    p = ModelParams(x_lo=0.2, x_hi=2.0, Nx=80)
    xg = x_grid(p)
    @test length(xg) == p.Nx
    @test xg[1] ≈ p.x_lo
    @test xg[end] ≈ p.x_hi
    @test all(diff(xg) .≈ (p.x_hi - p.x_lo) / (p.Nx - 1))
end

@testset "F_cdf and f_pdf" begin
    p = ModelParams(x_lo=0.2, x_hi=2.0)
    @test F_cdf(p.x_lo, p) ≈ 0.0
    @test F_cdf(p.x_hi, p) ≈ 1.0
    @test F_cdf((p.x_lo + p.x_hi) / 2, p) ≈ 0.5
    @test F_cdf(p.x_lo - 1.0, p) == 0.0
    @test F_cdf(p.x_hi + 1.0, p) == 1.0

    @test f_pdf(p) ≈ 1.0 / (p.x_hi - p.x_lo)
end

@testset "Continuation value" begin
    p = ModelParams()
    VS = 10.0
    VM_next = 15.0
    W = W_cont(VM_next, VS, p)
    @test W ≈ p.β * ((1.0 - p.λ) * VM_next + p.λ * VS)
end

@testset "find_x_cut" begin
    p = ModelParams()
    xg = x_grid(p)
    # VM0(x) = x so target recovers directly
    VM0_vals = collect(xg)
    target = 0.5 * (p.x_lo + p.x_hi)
    W = target * p.β - p.δ_C
    @test find_x_cut(W, VM0_vals, xg, p) ≈ target atol = 1e-6

    # Target below grid min: returns x_lo
    VM0_high = fill(100, length(xg))
    @test find_x_cut(-1000, VM0_high, xg, p) == p.x_lo

    # Target above grid max: returns x_hi
    VM0_low = fill(-100, length(xg))
    @test find_x_cut(1000, VM0_low, xg, p) == p.x_hi
end

@testset "compute_pc" begin
    p = ModelParams(α_prime=0.15)
    # for reservation value = x_lo cheater always leaves if found someone
    @test compute_pc(p.x_lo, p) ≈ p.α_prime
    # Nheater always leaves if found someone
    @test compute_pc(p.x_hi, p) ≈ 0.0
end

@testset "find_x_star" begin
    p = ModelParams()
    xg = x_grid(p)
    # VM0(x) = x, so x* ≈ VS
    VM0_vals = collect(xg)
    VS = 0.5 * (p.x_lo + p.x_hi)
    @test find_x_star(VM0_vals, VS, xg, p) ≈ VS atol = 1e-6

    # VS below grid
    @test find_x_star(VM0_vals, p.x_lo - 1.0, xg, p) == p.x_lo
    # VS above grid
    @test find_x_star(VM0_vals, p.x_hi + 1.0, xg, p) == p.x_hi
end

@testset "integral_VM0_above" begin
    p = ModelParams(x_lo=0.2, x_hi=2.0, Nx=200)
    xg = x_grid(p)
    density = f_pdf(p)

    # constant integrand
    c = 3.5
    VM0_const = fill(c, length(xg))
    a = 0.8
    I_const = integral_VM0_above(a, VM0_const, xg, p)
    @test I_const ≈ c * (p.x_hi - a) / (p.x_hi - p.x_lo) atol = 1e-6

    # x_cut at x_hi so integral is zero
    I_zero = integral_VM0_above(p.x_hi, VM0_lin, xg, p)
    @test I_zero ≈ 0.0 atol = 1e-10

    # x_cut below x_lo integral equals full-support expectation of the uniform distribution
    I_below = integral_VM0_above(p.x_lo - 1.0, VM0_lin, xg, p)
    @test I_below ≈ (p.x_lo + p.x_hi) / 2 atol = 1e-4

    # larger x_cut should yield smaller integral
    @test integral_VM0_above(0.5, VM0_lin, xg, p) > integral_VM0_above(1.5, VM0_lin, xg, p)
end
