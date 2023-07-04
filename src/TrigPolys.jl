module TrigPolys

using FFTW

export TrigPoly, random_trig_poly, evaluate, evaluateT, interpolatev, interpolate
export pad_to, truncate

"""
    struct TrigPoly{T<:AbstractFloat, VT<:AbstractVector{T}}
        a0::T    # Constant coefficient
        ac::VT   # cos coefficients
        as::VT   # sin coefficients
    end

Represents a *Hermitian trigonometric polynomial* by its coefficients.
The vectors `ac` and `as` should have the same length, that we call `n` in this
docstring.
This represent the following function
```julia
R(ω) = a0 + sum(ac[k] * cos(k * ω) + as[k] * sin(k * ω) for k in 1:n)
```
which is a polynomial in the variable `x = cos(ω)`.
"""
struct TrigPoly{T<:AbstractFloat, VT<:AbstractVector{T}}
    a0::T    # Constant coefficient
    ac::VT   # cos coefficients
    as::VT   # sin coefficients
    function TrigPoly(a0::S, ac::VS, as::VS) where {S<:AbstractFloat, VS<:AbstractVector{S}}
        @assert length(ac) == length(as) "sin and cos coefficients must have same length!"
        new{S, VS}(a0, ac, as)
    end
end

TrigPoly(a0::S, ac::VT, as::VU) where {
    S<:Number, T<:Number, U<:Number, VT<:AbstractVector{T}, VU<:AbstractVector{U}
} = TrigPoly(float(a0), float(ac), float(as))

TrigPoly(a::AbstractFloat) = TrigPoly(a, typeof(a)[], typeof(a)[])


"""
    TrigPoly(a::Number)

Converts a number `a` into a TrigPoly with a constant coefficient.
`a` will be converted to a float.
"""
TrigPoly(a::Number) = TrigPoly(float(a))

"""
    a0(p::TrigPoly)

Equivalent to `p.a0`, used for incorporation into automatic differentiation libraries.
"""
a0(p::TrigPoly) = p.a0

"""
    ac(p::TrigPoly)

Equivalent to `p.ac`, used for incorporation into automatic differentiation libraries.
"""
ac(p::TrigPoly) = p.ac

"""
    as(p::TrigPoly)

Equivalent to `p.as`, used for incorporation into automatic differentiation libraries.
"""
as(p::TrigPoly) = p.as

function Base.getproperty(p::TrigPoly, v::Symbol)
    if v == :n
        return length(p.ac)
    else
        return getfield(p, v)
    end
end

"""
    n(p::TrigPoly)

Equivalent to `p.n`. Returns the length of cosine/sine coefficients,
which is `length(p.ac) == length(p.as)`.
"""
n(p::TrigPoly) = p.n


"""
    degree(p::TrigPoly)

Returns the number of coefficients of `p`, which is `2*p.n + 1`.
"""
degree(p::TrigPoly) = 2 * p.n + 1

function get_mn(u::AbstractVector)
    m = length(u)
    @assert m % 2 == 1 "Only odd-length vectors supported!"
    n = (m - 1)÷2
    (m,n)
end

"""
    TrigPoly(u::AbstractVector{T}) where {T<:AbstractFloat}

Constructs a `TrigPoly` given a vector of coefficients of length `2n+1`. The
first entry of `u` is `a0`, the next `n` entries are `ac` and the last `n`
entries are `as`.
"""
function TrigPoly(u::AbstractVector{T}) where {T<:AbstractFloat}
    m, n = get_mn(u)
    TrigPoly(u[1], u[2:n+1], u[n+2:m])
end

"""
    Base.vec(p::TrigPoly)

Converts a TrigPoly into a vector, inverse operation of
[`TrigPolys.TrigPoly`](@ref TrigPolys.TrigPoly(u::AbstractVector{T}) where {T<:AbstractFloat})
"""
Base.vec(p::TrigPoly) = [a0(p); ac(p); as(p)]

"""
    Base.:(==)(p1::TrigPoly, p2::TrigPoly)

Returns true if the coefficients of `p1` and `p2` are all equal.
"""
function Base.:(==)(p1::TrigPoly, p2::TrigPoly)
    p1.n == p2.n && a0(p1) == a0(p2) && ac(p1) == ac(p2) && as(p1) == as(p2)
end

"""
    Base.:+(p1::TrigPoly, p2::TrigPoly)

Adds `p1` to `p2`. The degree of the resulting polynomial will be the maximum
degree of `p1` and `p2`: `(p1 + p2).n == max(p1.n, p2.n)`
"""
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

"""
    Base.:*(p1::TrigPoly, p2::TrigPoly)

Multiplies `p1` and `p2`. The degree of the resulting polynomial will be the sum of
degrees of `p1` and `p2`: `(p1 * p2).n == p1.n + p2.n`
"""
function Base.:*(p1::TrigPoly, p2::TrigPoly)
    n = p1.n + p2.n
    interpolate(evaluate(pad_to(p1, n)) .* evaluate(pad_to(p2, n)))
end
Base.:*(a::Number, p::TrigPoly)     = *(promote(p, a)...)
Base.:*(p::TrigPoly, a::Number)     = *(promote(p, a)...)

"""
    Base.:/(p::TrigPoly, a::Number)

Divides `p` by a scalar.
"""
Base.:/(p::TrigPoly, a::Number)     = p * (1/a)

Base.promote_rule(::Type{TrigPoly{S, VS}}, ::Type{T}) where {S<:AbstractFloat,VS<:AbstractVector{S},T<:Number} = TrigPoly{S, VS}
Base.convert(::Type{TrigPoly{S, VS}}, x::T) where {S<:AbstractFloat,VS<:AbstractVector{S},T<:Number} = TrigPoly(S(x))

"""
    pad_to(p::TrigPoly{T,VT}, n::Integer) where {T<:AbstractFloat, VT<:AbstractVector{T}}

Increase the length of `p.ac` and `p.as` to `n`, padding with zeros if necessary.
`n` cannot be smaller than `p.n`.
"""
function pad_to(p::TrigPoly{T,VT}, n::Integer) where {T<:AbstractFloat, VT<:AbstractVector{T}}
    @assert n >= p.n "Cannot pad to smaller degree!"
    TrigPoly(a0(p), VT([ac(p); zeros(n-p.n)]), VT([as(p); zeros(n-p.n)]))
end

function pad_to(v, n::Integer)
    xn = length(v)
    xs = (xn - 1) ÷ 2
    [v[1:xs+1]; zeros(n - xs); v[xs+2:xn]; zeros(n - xs)]
end

"""
    pad_by(p::TrigPoly, n::Integer)

Increase the length of `p.ac` and `p.as` by `n`, padding with `n` zeros.
"""
pad_by(p::TrigPoly, n::Integer) = pad_to(p, p.n + n)

"""
    truncate(p::TrigPoly, n::Integer)

Reduce the length of `p.ac` and `p.as` to `n` by deleting coefficients.
`n` cannot be larger than `p.n`.
"""
function truncate(p::TrigPoly, n::Integer)
    @assert n <= p.n "Cannot truncate to higher degree!"
    TrigPoly(a0(p), ac(p)[1:n], as(p)[1:n])
end


"""
    (p::TrigPoly)(x::Number)

Evaluate `p(x)` at a single point. See [`TrigPolys.evaluate`](@ref) for faster
evaluation for special points on the unit circle.
"""
function (p::TrigPoly)(x::Number)
    a0(p) + sum([ac(p)[k]*cos(k*x) + as(p)[k]*sin(k*x) for k in 1:p.n])
end

"""
    random_trig_poly(n)

Generate a `TrigPoly` with random coefficients.
"""
function random_trig_poly(n)
    TrigPoly(randn(), randn(n), randn(n))
end

"""
    basis(n, x)

Generate the basis for `TrigPoly` expressed as a vector, so that
`basis(p.n, x)'*vec(p) = p(x)`.
"""
basis(n, x) = vcat([1], [cos(k*x) for k in 1:n], [sin(k*x) for k in 1:n])


"""
    evaluate(u::AbstractVector)

Evaluates an array representation of `p::TrigPoly`, returned by `vec(p)`.
"""
function evaluate(u::AbstractVector)
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

"""
    evaluate(p::TrigPoly)

Evaluates `p` on `degree(p)` uniformly-spaced points on the circle. See
[here](https://en.wikipedia.org/wiki/Discrete_Fourier_transform#Trigonometric_interpolation_polynomial)
for more details.
"""
evaluate(p::TrigPoly) = evaluate(vec(p))

"""
    evaluateT(u::AbstractVector)

Adjoint of evaluate operator on an array representation of `p::TrigPoly`
, returned by `vec(p)`.
"""
function evaluateT(u::AbstractVector)
    m, n = get_mn(u)
    x = fft(u)
    [real(x[1:n+1]); reverse(imag(x[n+2:m]))]
end

"""
    evaluateT(p::TrigPoly)

Adjoint of evaluate operator. Returns a vector.
"""
evaluateT(p::TrigPoly) = evaluateT(vec(p))

"""
    interpolatev(u::Vector)

Inverse of [`TrigPolys.evaluate(p::AbstractVector)`](@ref). Returns a vector.
"""
function interpolatev(u::Vector)
    m, n = get_mn(u)
    out = evaluateT(u::Vector)
    for i = 1:m
        out[i] /= n+.5
    end
    out[1] /= 2
    out
end

"""
    interpolate(u::Vector)

Inverse of [`TrigPolys.evaluate(p::TrigPoly)`](@ref). Returns a `TrigPoly`.
"""
interpolate(u::Vector) = TrigPoly(interpolatev(u))

end
