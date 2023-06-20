# This file is k part of AffineMaps.jl, licensed under the MIT License (MIT).


"""
    struct AffineStep

Represents a multiplicative or additive step in an affine transformation
or it's inverse.

Multiplicative step are left-applied, additive steps are right-applied:

```julia
AffineStep(*, k)(x) ≈ k * x
AffineStep(+, k)(x) ≈ x + k
AffineStep(\\, k)(x) ≈ k \\ x
AffineStep(-, k)(x) ≈ x - k
```

`AffineStep(op, k)` supports `InverseFunctions.inverse` and
`ChangesOfVariables.with_logabsdet_jacobian`.
"""
struct AffineStep{F<:Union{typeof(*),typeof(\),typeof(+),typeof(-)},T} <: Function
    k::T
end

export AffineStep

AffineStep(f::F, k::T) where {F,T} = AffineStep{F,T}(k)

_fsym_char(::Type{typeof(*)}) = '*'
_fsym_char(::Type{typeof(\)}) = '\\'
_fsym_char(::Type{typeof(+)}) = '+'
_fsym_char(::Type{typeof(-)}) = '-'


function Base.show(io::IO, m::MIME"text/plain", f::AffineStep{F}) where F
    print(io, nameof(typeof(f)), "(", _fsym_char(F), ", ")
    show(io, m, f.k)
    print(io, ")")
end

function Base.show(io::IO, f::AffineStep{F}) where F
    print(io, nameof(typeof(f)), "(", _fsym_char(F), ", ")
    show(io, f.k)
    print(io, ")")
end


(f::AffineStep{typeof(*)})(x) = f.k * x
(f::AffineStep{typeof(\)})(x) = f.k \ x
(f::AffineStep{typeof(+)})(x) = x + f.k
(f::AffineStep{typeof(-)})(x) = x - f.k


InverseFunctions.inverse(f::AffineStep{typeof(+)}) = AffineStep(-, f.k)
InverseFunctions.inverse(f::AffineStep{typeof(-)}) = AffineStep(+, f.k)
InverseFunctions.inverse(f::AffineStep{typeof(*)}) = AffineStep(\, f.k)
InverseFunctions.inverse(f::AffineStep{typeof(\)}) = AffineStep(*, f.k)


# Julia v1.8 supports logabsdet(::Number), but older versions don't:
_logabsdet(x::Number) = log(abs(x))
_logabsdet(x::AbstractMatrix) = first(LinearAlgebra.logabsdet(x))

_type_ndof(::Type{<:Real}) = 1
_type_ndof(::Type{<:Complex}) = 2

_mul_ladj(k, x) = _logabsdet(k) * length(eachindex(x)) / length(axes(k,1)) * _type_ndof(eltype(x))

_realtype(::Type{T}) where {T<:Real} = T
_realtype(::Type{Complex{T}}) where {T<:Real} = T

const _RCNumber = Union{Real,Complex}

ChangesOfVariables.with_logabsdet_jacobian(f::AffineStep{typeof(+)}, x) = f(x), zero(_realtype(eltype(x)))
ChangesOfVariables.with_logabsdet_jacobian(f::AffineStep{typeof(-)}, x) = f(x), zero(_realtype(eltype(x)))

ChangesOfVariables.with_logabsdet_jacobian(f::AffineStep{typeof(*),<:Real}, x::Union{_RCNumber,Array{<:_RCNumber}}) = f(x), _mul_ladj(f.k, x)
ChangesOfVariables.with_logabsdet_jacobian(f::AffineStep{typeof(*),<:Complex}, x::Union{Complex,Array{<:Complex}}) = f(x), _mul_ladj(f.k, x)
ChangesOfVariables.with_logabsdet_jacobian(f::AffineStep{typeof(*),<:AbstractMatrix{<:Real}}, x::Union{AbstractVector{<:_RCNumber},AbstractMatrix{<:_RCNumber}}) = f(x), _mul_ladj(f.k, x)
ChangesOfVariables.with_logabsdet_jacobian(f::AffineStep{typeof(*),<:AbstractMatrix{<:Complex}}, x::Union{AbstractVector{<:Complex},AbstractMatrix{<:Complex}}) = f(x), _mul_ladj(f.k, x)

ChangesOfVariables.with_logabsdet_jacobian(f::AffineStep{typeof(\),<:Real}, x::Union{_RCNumber,Array{<:_RCNumber}}) = f(x), - _mul_ladj(f.k, x)
ChangesOfVariables.with_logabsdet_jacobian(f::AffineStep{typeof(\),<:Complex}, x::Union{Complex,Array{<:Complex}}) = f(x), - _mul_ladj(f.k, x)
ChangesOfVariables.with_logabsdet_jacobian(f::AffineStep{typeof(\),<:AbstractMatrix{<:Real}}, x::Union{AbstractVector{<:_RCNumber},AbstractMatrix{<:_RCNumber}}) = f(x), - _mul_ladj(f.k, x)
ChangesOfVariables.with_logabsdet_jacobian(f::AffineStep{typeof(\),<:AbstractMatrix{<:Complex}}, x::Union{AbstractVector{<:Complex},AbstractMatrix{<:Complex}}) = f(x), - _mul_ladj(f.k, x)
