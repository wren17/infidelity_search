using MarriageModel
using Test

@testset "GammaSpec" begin
    lin = LinearGamma()
    con = ConcaveGamma()
    @test gamma_val(lin, 0) == 0.0
    @test gamma_val(con, 0) == 0.0

    @test gamma_val(lin, 20) ≈ 1.6
    @test gamma_val(con, 100) ≈ 1.5 atol=0.1 # approaches γ_max

    #check concavity
    Δ1 = gamma_val(con, 2) - gamma_val(con, 1)
    Δ2 = gamma_val(con, 11) - gamma_val(con, 10)
    @test Δ1 > Δ2

end

@testset "DefaulParameters" begin
    p = ModelParams()
    @test p.δ_F <= p.δ_A <= p.δ_C
    @test p.α >= p.α_prime
    @test p.κ >= p.c
end