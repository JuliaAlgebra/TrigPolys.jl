# TrigPolys

[TrigPolys.jl](https://github.com/yuanchenyang/TrigPolys.jl) is a package for
fast manipulation of trigonometric polynomials.


A trignometric polynomial is defined on $$x \in [0,2\pi)$$ by

```math
p(x) = a_0 + \sum_{k=1}^n a_k \cos(kx) + a_{-k} \sin(kx)
```

The polynomial $$p(x)$$ can be represented either by $$2n+1$$ coefficients
$$a_k$$ or by evaluations at $$2n+1$$ distinct points in the interval
$$[0,2\pi)$$. Each representation is useful in different situations: the
coefficient representation is useful for truncating or increasing the degree of
a polynomial whereas the evaluation representation is useful for adding and
multiplying polynomials.

This package provides the functions `evaluate` and `interpolate` to convert
efficiently between these two representations. These operations are implemented
via the Fast Fourier Transform (FFT) provided by the
[FFTW.jl](https://github.com/JuliaMath/FFTW.jl) library.

For example, multiplying two trigonometric polynomials is implemented with the
following code:

```julia
function Base.:*(p1::TrigPoly, p2::TrigPoly)
    n = p1.n + p2.n
    interpolate(evaluate(pad_to(p1, n)) .* evaluate(pad_to(p2, n)))
end
```

## Construction and Convertion
```@docs
TrigPolys.TrigPoly
TrigPolys.TrigPoly(a::Number)
TrigPolys.TrigPoly(u::AbstractVector{T}) where {T<:AbstractFloat}
TrigPolys.a0
TrigPolys.ac
TrigPolys.as
TrigPolys.n
TrigPolys.degree
```

## Methods

```@docs
Base.vec(p::TrigPoly)
Base.:(==)(p1::TrigPoly, p2::TrigPoly)
Base.:+(p1::TrigPoly, p2::TrigPoly)
Base.:*(p1::TrigPoly, p2::TrigPoly)
Base.:/(p::TrigPoly, a::Number)
```

## Functions
```@docs
TrigPolys.pad_to
TrigPolys.pad_by
TrigPolys.random_trig_poly
TrigPolys.basis
TrigPolys.evaluate
TrigPolys.evaluateT
TrigPolys.interpolatev
TrigPolys.interpolate
```
