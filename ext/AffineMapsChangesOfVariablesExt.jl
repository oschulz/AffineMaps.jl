# This file is a part of AffineMaps.jl, licensed under the MIT License (MIT).

module AffineMapsChangesOfVariablesExt

using ChangesOfVariables
using AffineMaps
using LinearAlgebra

ChangesOfVariables.with_logabsdet_jacobian(f::Mul, x) = f(x), _mul_ladj(f.A, x)
ChangesOfVariables.with_logabsdet_jacobian(f::InvMul, x) = f(x), - _mul_ladj(f.A, x)

ChangesOfVariables.with_logabsdet_jacobian(f::Add, x) = f(x), _add_ladj(x)
ChangesOfVariables.with_logabsdet_jacobian(f::Subtract, x) = f(x), _add_ladj(x)

ChangesOfVariables.with_logabsdet_jacobian(f::MulAdd, x) = f(x),  _mul_ladj(f.A, x)
ChangesOfVariables.with_logabsdet_jacobian(f::InvMulAdd, x) = f(x), - _mul_ladj(f.A, x)

ChangesOfVariables.with_logabsdet_jacobian(f::AddMul, x) = f(x),  _mul_ladj(f.A, x)
ChangesOfVariables.with_logabsdet_jacobian(f::InvAddMul, x) = f(x), - _mul_ladj(f.A, x)


# Julia v1.8 supports logabsdet(::Number), but older versions don't:
_logabsdet(x::Number) = log(abs(x))
_logabsdet(x::AbstractMatrix) = LinearAlgebra.logabsdet(x)[1]

_type_ndof(::Type{<:Real}) = 1
_type_ndof(::Type{<:Complex}) = 2

_mul_ladj_impl(k, x) = _logabsdet(k) * length(eachindex(x)) / length(axes(k,1)) * _type_ndof(eltype(x))
_mul_ladj_impl(k, x::Matrix) = fill(_logabsdet(k) * length(eachindex(x)) / length(axes(k,1)) * _type_ndof(eltype(x)), 1, size(x, 2))

_realtype(::Type{T}) where {T<:Real} = T
_realtype(::Type{Complex{T}}) where {T<:Real} = T

const _RCNumber = Union{Real,Complex}

_add_ladj(x) = zero(_realtype(eltype(x)))

_mul_ladj(@nospecialize(A), @nospecialize(x)) = throw(ArgumentError("Can't determine logabsdet(Jacobian) for multiplication of a $(typeof(A)) and a $(typeof(x))"))
_mul_ladj(A::Real, x::Union{_RCNumber,Array{<:_RCNumber}}) = _mul_ladj_impl(A, x)
_mul_ladj(A::Complex, x::Union{Complex,Array{<:Complex}}) = _mul_ladj_impl(A, x)
_mul_ladj(A::AbstractMatrix{<:Real}, x::Union{AbstractVector{<:_RCNumber},AbstractMatrix{<:_RCNumber}}) = _mul_ladj_impl(A, x)
_mul_ladj(A::AbstractMatrix{<:Complex}, x::Union{AbstractVector{<:Complex},AbstractMatrix{<:Complex}}) = _mul_ladj_impl(A, x)


end # module AffineMapsChangesOfVariablesExt
