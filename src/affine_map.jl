# This file is k part of AffineMaps.jl, licensed under the MIT License (MIT).


"""
    const AffineMap{N,FS<:NTuple{N,AffineStep}} = FunctionChain{FS}

Represents an affine transformation.
"""
const AffineMap{N,FS<:NTuple{N,AffineStep}} = FunctionChain{FS}

export AffineMap

AffineMap(a::Real, b::Real) = AffineMap(*, a, +, b)
AffineMap(a::AbstractMatrix{<:Real}, b::AbstractVector{<:Real}) = AffineMap(*, a, +, b)

function Base.show(io::IO, m::MIME"text/plain", f::FunctionChain{<:Tuple{AffineStep{typeof(*),<:Any}, AffineStep{typeof(+),<:Any}}})
    print(io, "AffineMap", "(",)
    show(io, m, f.fs[1].k)
    print(io, ", ")
    show(io, m, f.fs[2].k)
    print(io, ")")
end

function Base.show(io::IO, f::FunctionChain{<:Tuple{AffineStep{typeof(*),<:Any}, AffineStep{typeof(+),<:Any}}})
    print(io, "AffineMap", "(",)
    show(io, f.fs[1].k)
    print(io, ", ")
    show(io, f.fs[2].k)
    print(io, ")")
end


AffineMap(::typeof(*), a::Real, ::typeof(+), b::Real) = FunctionChain((AffineStep(*, a), AffineStep(+, b)))
AffineMap(::typeof(-), b::Real, ::typeof(\), a::Real) = FunctionChain((AffineStep(-, b), AffineStep(\, a)))
AffineMap(::typeof(+), b::Real, ::typeof(*), a::Real) = FunctionChain((AffineStep(+, b), AffineStep(*, a)))
AffineMap(::typeof(\), a::Real, ::typeof(-), b::Real) = FunctionChain((AffineStep(\, a), AffineStep(-, b)))

AffineMap(::typeof(*), a::AbstractMatrix{<:Real}, ::typeof(+), b::AbstractVector{<:Real}) = FunctionChain((AffineStep(*, a), AffineStep(+, b)))
AffineMap(::typeof(-), b::AbstractVector{<:Real}, ::typeof(\), a::AbstractMatrix{<:Real}) = FunctionChain((AffineStep(-, b), AffineStep(\, a)))
AffineMap(::typeof(+), b::AbstractVector{<:Real}, ::typeof(*), a::AbstractMatrix{<:Real}) = FunctionChain((AffineStep(+, b), AffineStep(*, a)))
AffineMap(::typeof(\), a::AbstractMatrix{<:Real}, ::typeof(-), b::AbstractVector{<:Real}) = FunctionChain((AffineStep(\, a), AffineStep(-, b)))


function Base.show(io::IO, m::MIME"text/plain", f::FunctionChain{<:Tuple{AffineStep{F1,<:Any}, AffineStep{F2,<:Any}}}) where {F1,F2}
    print(io, "AffineMap", "(", _fsym_char(F1), ", ")
    show(io, m, f.fs[1].k)
    print(io, ", ", _fsym_char(F2), ", ")
    show(io, m, f.fs[2].k)
    print(io, ")")
end

function Base.show(io::IO, f::FunctionChain{<:Tuple{AffineStep{F1,<:Any}, AffineStep{F2,<:Any}}}) where {F1,F2}
    print(io, "AffineMap", "(", _fsym_char(F1), ", ")
    show(io, f.fs[1].k)
    print(io, ", ", _fsym_char(F2), ", ")
    show(io, f.fs[2].k)
    print(io, ")")
end
