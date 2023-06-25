# This file is k part of AffineMaps.jl, licensed under the MIT License (MIT).

"""
    abstract type AbstractAffineMap <: Function

Abstract type for affine maps.

Affine map `f::AbstractAffineMap` act like `f(x) == A * x .+ b` or
`f(x) == A * (x .+ b)` or their inverses, depending on the type of
the map `f`.

`A` may be a `Number`,  `AbstractMatrix{<:Number}` or any other
multiplicative linear operator in general (that supports at least `size(A)`
and `eltype(A)`). The packages
[LinearMaps](https://github.com/JuliaLinearAlgebra/LinearMaps.jl),
[LinearOperators](https://github.com/JuliaSmoothOptimizers/LinearOperators.jl)
and [SciMLOperators](https://github.com/SciML/SciMLOperators.jl)
provide such operators, for example.

`b` must have a shape that supports broadcasted addition and subtraction with
`x`.

Subtypes of `AbstractAffineMap` should implement/support the APIs of

* [InverseFunctions.jl](https://github.com/JuliaMath/InverseFunctions.jl)
* [ChangesOfVariables.jl](https://github.com/JuliaMath/ChangesOfVariables.jl)
* [Functors.jl](https://github.com/FluxML/Functors.jl)
"""
abstract type AbstractAffineMap <: Function end


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

Base.:(==)(f::Mul, g::Mul) = f.A == g.A
Base.isapprox(f::Mul, g::Mul; kwargs...) = isapprox(f.A, g.A; kwargs...)

"""
    struct Add

`f = Add(b)` has the behavior `f(x) == x .+ f.b`.

See [`AbstractAffineMap`](@ref) for more information.
"""
struct Add{V} <: AbstractAffineMap
    b::V
end
export Add

(f::Add)(x) = _bc_add(x, f.b)

Base.:(==)(f::Add, g::Add) = f.b == g.b
Base.isapprox(f::Add, g::Add; kwargs...) = isapprox(f.b, g.b; kwargs...)


"""
    struct MulAdd

`f = MulAdd(A, b)` has the behavior `f(x) == f.A * x .+ f.b`.

See [`AbstractAffineMap`](@ref) for more information.
"""
struct MulAdd{M,V} <: AbstractAffineMap
    A::M
    b::V
end
export MulAdd

(f::MulAdd)(x) = _bc_muladd(f.A, x, f.b)

Base.:(==)(f::MulAdd, g::MulAdd) = f.A == g.A && f.b == g.b
Base.isapprox(f::MulAdd, g::MulAdd; kwargs...) = isapprox(f.A, g.A; kwargs...) && isapprox(f.b, g.b; kwargs...)


"""
    struct AddMul

`f = AddMul(A, b)` has the behavior `f(x) == f.A * (x .+ f.b)`.

See [`AbstractAffineMap`](@ref) for more information.
"""
struct AddMul{V,M} <: AbstractAffineMap
    b::V
    A::M
end
export AddMul

(f::AddMul)(x) = f.A * _bc_add(x, f.b)

Base.:(==)(f::AddMul, g::AddMul) = f.A == g.A && f.b == g.b
Base.isapprox(f::AddMul, g::AddMul; kwargs...) = isapprox(f.A, g.A; kwargs...) && isapprox(f.b, g.b; kwargs...)


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

Base.:(==)(f::InvMul, g::InvMul) = f.A == g.A
Base.isapprox(f::InvMul, g::InvMul; kwargs...) = isapprox(f.A, g.A; kwargs...)


"""
    struct Subtract

`f = Subtract(b)` has the behavior `f(x) == x .- f.b`. It is the inverse of
`Add(b)`.

See [`AbstractAffineMap`](@ref) for more information.
"""
struct Subtract{V} <: AbstractAffineMap
    b::V
end
export Subtract

(f::Subtract)(x) = _bc_sub(x, f.b)

Base.:(==)(f::Subtract, g::Subtract) = f.b == g.b
Base.isapprox(f::Subtract, g::Subtract; kwargs...) = isapprox(f.b, g.b; kwargs...)


"""
    struct InvMulAdd

`f = InvMulAdd(A, b)` has the behavior `f(x) == f.A \\ (x .- f.b)`. It is the
inverse of `MulAdd(A, b)`.

See [`AbstractAffineMap`](@ref) for more information.
"""
struct InvMulAdd{M,V} <: AbstractAffineMap
    A::M
    b::V
end
export InvMulAdd

(f::InvMulAdd)(x) = f.A \ _bc_sub(x, f.b)

Base.:(==)(f::InvMulAdd, g::InvMulAdd) = f.A == g.A && f.b == g.b
Base.isapprox(f::InvMulAdd, g::InvMulAdd; kwargs...) = isapprox(f.A, g.A; kwargs...) && isapprox(f.b, g.b; kwargs...)


""" muladd(f.A, x, f.b)
    struct InvAddMul

`f = InvAddMul(A, b)` has the behavior `f(x) == (f.A \\ x) .- f.b`. It is the
inverse of `AddMul(A, b)`.

See [`AbstractAffineMap`](@ref) for more information.
"""
struct InvAddMul{V,M} <: AbstractAffineMap
    b::V
    A::M
end
export InvAddMul

(f::InvAddMul)(x) = _bc_sub(f.A \ x, f.b)

Base.:(==)(f::InvAddMul, g::InvAddMul) = f.A == g.A && f.b == g.b
Base.isapprox(f::InvAddMul, g::InvAddMul; kwargs...) = isapprox(f.A, g.A; kwargs...) && isapprox(f.b, g.b; kwargs...)


_bc_add(a::Number, b::Number) = a + b
_bc_add(a::AbstractArray{T,N}, b::AbstractArray{U,N}) where {T,U,N} = a + b
_bc_add(a, b) = a .+ b

_bc_sub(a::Number, b::Number) = a - b
_bc_sub(a::AbstractArray{T,N}, b::AbstractArray{U,N}) where {T,U,N} = a - b
_bc_sub(a, b) = a .- b

_bc_muladd(a, b, c) = _bc_add(a * b, c)
_bc_muladd(a::Number, b::Number, c::Number) = muladd(a, b, c)
_bc_muladd(a::Number, b::AbstractVector, c::AbstractVector) = muladd(a, b, c)
_bc_muladd(a::Number, b::AbstractMatrix, c::AbstractMatrix) = muladd(a, b, c)
_bc_muladd(a::AbstractMatrix, b::Number, c::AbstractMatrix) = muladd(a, b, c)
_bc_muladd(a::AbstractMatrix, b::AbstractVector, c::Number) = muladd(a, b, c)
_bc_muladd(a::AbstractMatrix, b::AbstractVector, c::AbstractVector) = muladd(a, b, c)
_bc_muladd(a::AbstractMatrix, b::AbstractMatrix, c::Number) = muladd(a, b, c)
_bc_muladd(a::AbstractMatrix, b::AbstractMatrix, c::AbstractVector) = muladd(a, b, c)
_bc_muladd(a::AbstractMatrix, b::AbstractMatrix, c::AbstractMatrix) = muladd(a, b, c)
