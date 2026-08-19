# This file is a part of AffineMaps.jl, licensed under the MIT License (MIT).

module AffineMapsChangesOfVariablesExt

using ChangesOfVariables
using AffineMaps
using LinearAlgebra

# The ladj computation carries the shape and element-type guards, so it
# runs before the primal (arguments evaluate before the call), turning
# invalid argument shapes into a uniform ArgumentError instead of a
# factor-type-dependent primal error or silent nonsense:
_primal_with_ladj(f, x, ladj) = f(x), ladj

ChangesOfVariables.with_logabsdet_jacobian(f::Mul, x) = _primal_with_ladj(f, x, _mul_ladj(f.A, x))
ChangesOfVariables.with_logabsdet_jacobian(f::InvMul, x) = _primal_with_ladj(f, x, - _mul_ladj(f.A, x))

ChangesOfVariables.with_logabsdet_jacobian(f::Add, x) = f(x), _add_ladj(x)
ChangesOfVariables.with_logabsdet_jacobian(f::Subtract, x) = f(x), _add_ladj(x)

ChangesOfVariables.with_logabsdet_jacobian(f::MulAdd, x) = _primal_with_ladj(f, x, _mul_ladj(f.A, x))
ChangesOfVariables.with_logabsdet_jacobian(f::InvMulAdd, x) = _primal_with_ladj(f, x, - _mul_ladj(f.A, x))

ChangesOfVariables.with_logabsdet_jacobian(f::AddMul, x) = _primal_with_ladj(f, x, _mul_ladj(f.A, x))
ChangesOfVariables.with_logabsdet_jacobian(f::InvAddMul, x) = _primal_with_ladj(f, x, - _mul_ladj(f.A, x))


_logabsdet(x) = first(LinearAlgebra.logabsdet(x))

_type_ndof(::Type{<:Real}) = 1
_type_ndof(::Type{<:Complex}) = 2

_mul_ladj_impl(A, x) = _logabsdet(A) * length(eachindex(x)) / size(A, 1) * _type_ndof(eltype(x))
_mul_ladj_impl(A, x::AbstractMatrix) = fill(_logabsdet(A) * size(x, 1) / size(A, 1) * _type_ndof(eltype(x)), 1, size(x, 2))

_realtype(::Type{T}) where {T<:Real} = T
_realtype(::Type{Complex{T}}) where {T<:Real} = T

const _RCNumber = Union{Real,Complex}

_add_ladj(x) = zero(_realtype(eltype(x)))
_add_ladj(x::AbstractMatrix) = zeros(_realtype(eltype(x)), 1, size(x, 2))

# Holy-trait dispatch on dimensionality, for factor types that carry a
# shape but no AbstractArray supertype:
struct _NDims{N} end
_NDims(x) = _NDims{ndims(x)}()

# The multiplicative factor of an affine map is linear by contract, so any
# matrix-shaped factor type that supports the LinearAlgebra interface
# (multiplication, ndims, size, eltype, logabsdet) has a constant
# log-abs-det Jacobian - linear operator types compose naturally, without
# dedicated glue code. Only matrix-shaped factors (ndims two) acting on
# vectors or column batches qualify, and real-valued factors act on real
# or complex spaces while complex-valued factors require complex arguments
# (a complex factor on a real space is not a change of variables within
# that space). Factor types without a logabsdet method fail with a
# MethodError naming the missing capability. The shape and element-type
# checks happen below the fully generic method, they must not narrow its
# signature (that would be ambiguous against the scalar-factor methods):
_mul_ladj(A, x) = _generic_mul_ladj(_NDims(A), A, x)

_generic_mul_ladj(::_NDims{2}, A, x::AbstractVecOrMat{<:_RCNumber}) = _mul_ladj_checked(eltype(A), eltype(x), A, x)
_generic_mul_ladj(::_NDims, @nospecialize(A), @nospecialize(x)) = throw(ArgumentError("Can't determine logabsdet(Jacobian) for multiplication of a $(typeof(A)) and a $(typeof(x))"))

_mul_ladj_checked(::Type{<:Real}, ::Type{<:_RCNumber}, A, x) = _mul_ladj_impl(A, x)
_mul_ladj_checked(::Type{<:Complex}, ::Type{<:Complex}, A, x) = _mul_ladj_impl(A, x)
_mul_ladj_checked(::Type, ::Type, @nospecialize(A), @nospecialize(x)) = throw(ArgumentError("Can't determine logabsdet(Jacobian) for multiplication of a $(typeof(A)) and a $(typeof(x))"))
_mul_ladj(A::UniformScaling, x) = _mul_ladj(A.λ, x)
_mul_ladj(A::Real, x::Union{_RCNumber,AbstractArray{<:_RCNumber}}) = _mul_ladj_impl(A, x)
_mul_ladj(A::Complex, x::Union{Complex,AbstractArray{<:Complex}}) = _mul_ladj_impl(A, x)
_mul_ladj(A::AbstractMatrix{<:Real}, x::Union{AbstractVector{<:_RCNumber},AbstractMatrix{<:_RCNumber}}) = _mul_ladj_impl(A, x)
_mul_ladj(A::AbstractMatrix{<:Complex}, x::Union{AbstractVector{<:Complex},AbstractMatrix{<:Complex}}) = _mul_ladj_impl(A, x)


end # module AffineMapsChangesOfVariablesExt
