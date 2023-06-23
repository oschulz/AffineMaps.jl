# This file is k part of AffineMaps.jl, licensed under the MIT License (MIT).

"""
    abstract type AbstractAffineMap

Abstract type for affine maps.

Affine map `f::AbstractAffineMap` act like `f(x) == A * x + b` or
`f(x) == A * (x + b)` or their inverses, depending on the type of
the map `f`.

`A` may be a `Number`,  `AbstractMatrix{<:Number}` or any other
multiplicative linear operator in general (that supports at least `size(A)`
and `eltype(A)`). The packages
[LinearMaps](https://github.com/JuliaLinearAlgebra/LinearMaps.jl),
[LinearOperators](https://github.com/JuliaSmoothOptimizers/LinearOperators.jl)
and [SciMLOperators](https://github.com/SciML/SciMLOperators.jl)
provide such operators, for example.

`b` must support addition and subtraction with `x`, so it may be
a `Number` or `AbstractVector{<:Number}`, depending on the type of `x`.

Subtypes of `AbstractAffineMap` should implement/support the APIs of

* [InverseFunctions.jl](https://github.com/JuliaMath/InverseFunctions.jl)
* [ChangesOfVariables.jl](https://github.com/JuliaMath/ChangesOfVariables.jl)
* [Functors.jl](https://github.com/FluxML/Functors.jl)
"""
abstract type AbstractAffineMap end


"""
    struct Mul

`f = Mul(A)` has the behavior `f(x) == f.A * x`.

See [`AbstractAffineMap`](@ref) for more information.
"""
struct Mul{M} <: AbstractAffineMap
    A::M
end
export Mul

(f::Mul)(x) = f.A * x


"""
    struct Add

`f = Add(b)` has the behavior `f(x) == x + f.b`.

See [`AbstractAffineMap`](@ref) for more information.
"""
struct Add{V} <: AbstractAffineMap
    b::V
end
export Add

(f::Add)(x) = x + f.b


"""
    struct MulAdd

`f = MulAdd(A, b)` has the behavior `f(x) == f.A * x + f.b`.

See [`AbstractAffineMap`](@ref) for more information.
"""
struct MulAdd{M,V} <: AbstractAffineMap
    A::M
    b::V
end
export MulAdd

(f::MulAdd)(x) = muladd(f.A, x, f.b)


"""
    struct AddMul

`f = AddMul(A, b)` has the behavior `f(x) == f.A * (x + f.b)`.

See [`AbstractAffineMap`](@ref) for more information.
"""
struct AddMul{V,M} <: AbstractAffineMap
    b::V
    A::M
end
export AddMul

(f::AddMul)(x) = f.A * (x + f.b)


"""
    struct InvMul

`f = InvMul(A)` has the behavior `f(x) == f.A \\ x`. It is the inverse of
`Mul(A)`.

See [`AbstractAffineMap`](@ref) for more information.
"""
struct InvMul{M} <: AbstractAffineMap
    A::M
end
export InvMul

(f::InvMul)(x) = f.A \ x


"""
    struct Subtract

`f = Subtract(b)` has the behavior `f(x) == x - f.b`. It is the inverse of
`Add(b)`.

See [`AbstractAffineMap`](@ref) for more information.
"""
struct Subtract{V} <: AbstractAffineMap
    b::V
end
export Subtract

(f::Subtract)(x) = x - f.b


"""
    struct InvMulAdd

`f = InvMulAdd(A, b)` has the behavior `f(x) == f.A \\ (x - f.b)`. It is the
inverse of `MulAdd(A, b)`.

See [`AbstractAffineMap`](@ref) for more information.
"""
struct InvMulAdd{M,V} <: AbstractAffineMap
    A::M
    b::V
end
export InvMulAdd

(f::InvMulAdd)(x) = f.A \ (x - f.b)


""" muladd(f.A, x, f.b)
    struct InvAddMul

`f = InvAddMul(A, b)` has the behavior `f(x) == (f.A \\ x) - f.b`. It is the
inverse of `AddMul(A, b)`.

See [`AbstractAffineMap`](@ref) for more information.
"""
struct InvAddMul{V,M} <: AbstractAffineMap
    b::V
    A::M
end
export InvAddMul

(f::InvAddMul)(x) = (f.A \ x) - f.b


InverseFunctions.inverse(f::Mul) = InvMul(f.A)
InverseFunctions.inverse(f::InvMul) = Mul(f.A)

InverseFunctions.inverse(f::Add) = Subtract(f.b)
InverseFunctions.inverse(f::Subtract) = Add(f.b)

InverseFunctions.inverse(f::MulAdd) = InvMulAdd(f.A, f.b)
InverseFunctions.inverse(f::InvMulAdd) = MulAdd(f.A, f.b)

InverseFunctions.inverse(f::AddMul) = InvAddMul(f.b, f.A)
InverseFunctions.inverse(f::InvAddMul) = AddMul(f.b, f.A)


ChangesOfVariables.with_logabsdet_jacobian(f::Mul, x) = f(x), _mul_ladj(f.A, x)
ChangesOfVariables.with_logabsdet_jacobian(f::InvMul, x) = f(x), - _mul_ladj(f.A, x)

ChangesOfVariables.with_logabsdet_jacobian(f::Add, x) = f(x), _add_ladj(x)
ChangesOfVariables.with_logabsdet_jacobian(f::Subtract, x) = f(x), _add_ladj(x)

ChangesOfVariables.with_logabsdet_jacobian(f::MulAdd, x) = f(x),  _mul_ladj(f.A, x)
ChangesOfVariables.with_logabsdet_jacobian(f::InvMulAdd, x) = f(x), - _mul_ladj(f.A, x)

ChangesOfVariables.with_logabsdet_jacobian(f::AddMul, x) = f(x),  _mul_ladj(f.A, x)
ChangesOfVariables.with_logabsdet_jacobian(f::InvAddMul, x) = f(x), -  _mul_ladj(f.A, x)


# Julia v1.8 supports logabsdet(::Number), but older versions don't:
_logabsdet(x::Number) = log(abs(x))
_logabsdet(x::AbstractMatrix) = LinearAlgebra.logabsdet(x)[1]

_type_ndof(::Type{<:Real}) = 1
_type_ndof(::Type{<:Complex}) = 2

_mul_ladj_impl(k, x) = _logabsdet(k) * length(eachindex(x)) / length(axes(k,1)) * _type_ndof(eltype(x))

_realtype(::Type{T}) where {T<:Real} = T
_realtype(::Type{Complex{T}}) where {T<:Real} = T

const _RCNumber = Union{Real,Complex}

_add_ladj(x) = zero(_realtype(eltype(x)))

_mul_ladj(@nospecialize(A), @nospecialize(x)) = throw(ArgumentError("Can't determine logabsdet(Jacobian) for multiplication a $(typeof(A)) and a $(typeof(x))"))
_mul_ladj(A::Real, x::Union{_RCNumber,Array{<:_RCNumber}}) = _mul_ladj_impl(A, x)
_mul_ladj(A::Complex, x::Union{Complex,Array{<:Complex}}) = _mul_ladj_impl(A, x)
_mul_ladj(A::AbstractMatrix{<:Real}, x::Union{AbstractVector{<:_RCNumber},AbstractMatrix{<:_RCNumber}}) = _mul_ladj_impl(A, x)
_mul_ladj(A::AbstractMatrix{<:Complex}, x::Union{AbstractVector{<:Complex},AbstractMatrix{<:Complex}}) = _mul_ladj_impl(A, x)


Functors.@functor Mul
Functors.@functor Add
Functors.@functor MulAdd
Functors.@functor AddMul
Functors.@functor InvMul
Functors.@functor Subtract
Functors.@functor InvMulAdd
Functors.@functor InvAddMul
