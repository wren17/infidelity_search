using Test
using MarriageModel

p = ModelParams()
sol = solve_model(p)


@testset "output shapes" begin
    Nx = length(sol.xg)
    @test size(sol.VM) == (p.T + 1, Nx)
    @test size(sol.eq_type_grid) == (p.T + 1, Nx)
    @test size(sol.s_star_grid) == (p.T + 1, Nx)
    @test length(sol.VM0_vals) == Nx
end

@testset "VM0 increasing in x" begin
    @test all(diff(sol.VM0_vals) .>= -1e-10)
end

@testset "x_star in grid range" begin
    @test sol.x_star >= p.x_lo
    @test sol.x_star <= p.x_hi
end


@testset "outer loop: initial guess VS updates" begin
    xg = x_grid(p)
    VS_init = p.b / (1.0 - p.β)
    VM0_init = [match_value(0, x, p) / (1.0 - p.β) for x in xg]
    VS_new = compute_VS(VM0_init, xg, VS_init, p)
    @test isfinite(VS_new)
    @test VS_new >= VS_init
end


