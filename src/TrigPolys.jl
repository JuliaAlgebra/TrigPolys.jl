module TrigPolys

using FFTW

export TrigPoly, random_trig_poly, evaluate, evaluateT, interpolatev, interpolate
export pad_to, truncate

# Represents a trigonometric polynomial by its coefficients
struct TrigPoly{T<:AbstractFloat, VT<:AbstractVector{T}}
    a0::T    # Constant coefficient
    ac::VT   # cos coefficients
    as::VT   # sin coefficients
    function TrigPoly(a0::S, ac::VS, as::VS) where {S <:AbstractFloat, VS<:AbstractVector{S}}
        @assert length(ac) == length(as) "sin and cos coefficients must have same length!"
        new{S, VS}(a0, ac, as)
    end
end

TrigPoly(a::AbstractFloat) = TrigPoly(a, typeof(a)[], typeof(a)[])
TrigPoly(a::Number) = TrigPoly(float(a))

a0(p::TrigPoly) = p.a0
ac(p::TrigPoly) = p.ac
as(p::TrigPoly) = p.as

degree(p::TrigPoly) = 2 * p.n + 1

function Base.getproperty(p::TrigPoly, v::Symbol)
    if v == :n
        return length(p.ac)
    else
        return getfield(p, v)
    end
end

function get_mn(u::AbstractArray)
    m = length(u)
    @assert m % 2 == 1 "Only odd-length vectors supported!"
    n = (m - 1)÷2
    (m,n)
end

function TrigPoly(u::AbstractArray)
    m, n = get_mn(u)
    TrigPoly(u[1], u[2:n+1], u[n+2:m])
end

Base.vec(p::TrigPoly) = [a0(p); ac(p); as(p)]

function Base.:(==)(p1::TrigPoly, p2::TrigPoly)
    p1.n == p2.n && a0(p1) == a0(p2) && ac(p1) == ac(p2) && as(p1) == as(p2)
end

function Base.:+(p1::TrigPoly, p2::TrigPoly)
    n = max(p1.n, p2.n)
    p1 = pad_to(p1, n)
    p2 = pad_to(p2, n)
    TrigPoly(a0(p1)+a0(p2), ac(p1)+ac(p2), as(p1)+as(p2))
end
Base.:+(p::TrigPoly, a::Number)     = +(promote(p, a)...)
Base.:+(a::Number, p::TrigPoly)     = p + a

Base.:-(p::TrigPoly)                = TrigPoly(-a0(p), -ac(p), -as(p))
Base.:-(p1::TrigPoly, p2::TrigPoly) = p1 + (-p2)
Base.:-(p::TrigPoly, a::Number)     = -(promote(p, a)...)
Base.:-(a::Number, p::TrigPoly)     = -(promote(p, a)...)

function Base.:*(p1::TrigPoly, p2::TrigPoly)
    n = p1.n + p2.n
    interpolate(evaluate(pad_to(p1, n)) .* evaluate(pad_to(p2, n)))
end
Base.:*(a::Number, p::TrigPoly)     = *(promote(p, a)...)
Base.:*(p::TrigPoly, a::Number)     = *(promote(p, a)...)

Base.:/(p::TrigPoly, a::Number)     = p * (1/a)

Base.promote_rule(::Type{TrigPoly{S, VS}}, ::Type{T}) where {S<:AbstractFloat,VS<:AbstractVector{S},T<:Number} = TrigPoly{S, VS}
Base.convert(::Type{TrigPoly{S, VS}}, x::T) where {S<:AbstractFloat,VS<:AbstractVector{S},T<:Number} = TrigPoly(x)

# Increase the maxdegree of a TrigPoly **to** n
function pad_to(p::TrigPoly, n::Integer)
    @assert n >= p.n "Cannot pad to smaller degree!"
    TrigPoly(a0(p), [ac(p); zeros(n-p.n)], [as(p); zeros(n-p.n)])
end

function pad_to(v, n::Integer)
    xn = length(v)
    xs = (xn - 1) ÷ 2
    [v[1:xs+1]; zeros(n - xs); v[xs+2:xn]; zeros(n - xs)]
end

# Increase the maxdegree of a TrigPoly **by** n
pad_by(p::TrigPoly, n::Integer) = pad_to(p, p.n + n)

# Reduce the maxdegree of a TrigPoly **to** n
function truncate(p::TrigPoly, n::Integer)
    @assert n <= p.n "Cannot truncate to higher degree!"
    TrigPoly(a0(p), ac(p)[1:n], as(p)[1:n])
end

function (p::TrigPoly)(x::Number)
    a0(p) + sum([ac(p)[k]*cos(k*x) + as(p)[k]*sin(k*x) for k in 1:p.n])
end

function random_trig_poly(n)
    TrigPoly(randn(), randn(n), randn(n))
end

basis(n, x) = vcat([1], [cos(k*x) for k in 1:n], [sin(k*x) for k in 1:n])

# https://en.wikipedia.org/wiki/Discrete_Fourier_transform#Trigonometric_interpolation_polynomial
function evaluate(u::AbstractArray)
    m, n = get_mn(u)
    out = similar(u, Complex{Float64})
    out[1] = u[1]
    for (i, (up, un)) in enumerate(zip(2:n+1, n+2:m))
       out[1+i]   = (u[up]+im*u[un])/2
       out[m+1-i] = (u[up]-im*u[un])/2
    end
    fft!(out)
    real(out)
    # Slower non-inplace version
    #   u0 = u[1]
    #   up = u[2:n+1]
    #   un = u[n+2:m]
    #   pos = (up-im*un)./2
    #   neg = (up+im*un)./2
    #   real(fft([u0; neg; reverse(pos)]))'
end

# Evaluate a trigonometric polynomial on uniformly-spaced points on the circle
evaluate(p::TrigPoly) = evaluate(vec(p))

# Adjoint of evaluate operator
function evaluateT(u::AbstractArray)
    m, n = get_mn(u)
    x = fft(u)
    [real(x[1:n+1]); reverse(imag(x[n+2:m]))]
end

evaluateT(p::TrigPoly) = evaluateT(vec(p))

# Inverse of evaluate operator
function interpolatev(u::Vector)
    m, n = get_mn(u)
    out = evaluateT(u::Vector)
    for i = 1:m
        out[i] /= n+.5
    end
    out[1] /= 2
    out
end

# Returns TrigPoly
interpolate(u::Vector) = TrigPoly(interpolatev(u))


end
