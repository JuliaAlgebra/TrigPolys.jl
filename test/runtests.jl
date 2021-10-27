using Test
using TrigPolys

@testset "Arithmetic" begin
    p1 = random_trig_poly(102)
    p2 = random_trig_poly(201)
    x = range(0; stop=2*pi, length=2000)

    @test isapprox((-p1).(x), -(p1.(x)))
    @test isapprox((p1 + p2).(x), p1.(x) .+ p2.(x))
    @test isapprox((p1 - p2).(x), p1.(x) .- p2.(x))
    @test isapprox((p1 * p2).(x), p1.(x) .* p2.(x))

    @test TrigPolys.a0(p1 + 3) == TrigPolys.a0(p1) + 3
    @test TrigPoly(pi) == TrigPoly(pi)
    @test p1 + p2 == p2 + p1
    @test p1 * p2 == p2 * p1
    @test p1 + pi == pi + p1
    @test p1 * pi == pi * p1
    @test p1 / pi == p1 * (1/pi)
end

@testset "Type Conversion and Promotion" begin
    F1 = Float32
    F2 = Float64
    VF1 = Vector{F1}
    VF2 = Vector{F2}
    @test typeof(TrigPoly(F2(4))+F1(3)) == TrigPoly{F2, VF2}
    @test typeof(convert(TrigPoly{F2, VF2}, F1(.3))) == TrigPoly{F2, VF2}
    @test typeof(convert(TrigPoly{F1, VF1}, F2(.3))) == TrigPoly{F1, VF1}
    @test typeof(TrigPoly(F1(2))+TrigPoly(F2(3))) == TrigPoly{F2, VF2}
    @test typeof(TrigPoly(F1(2))+F1(3)) == TrigPoly{F1, VF1}
end

@testset "Evaluate and interpolate" begin
    n = 1000
    p = random_trig_poly(n)
    s = p.n*2+1
    samples = [2*pi*(i-1)/s for i in 1:s]
    B = hcat([TrigPolys.basis(p.n, xi) for xi in samples]...)

    u = randn(s)
    @test isapprox(evaluate(u), B'*u)
    @test isapprox(evaluate(p), p.(samples))
    @test isapprox(evaluateT(vec(u)), B*u)
    @test isapprox(interpolatev(evaluate(u)), u)
end
