# This file is a part of AffineMaps.jl, licensed under the MIT License (MIT).

module AffineMapsFunctorsExt

using Functors
using AffineMaps

@static if !isdefined(Base, :pkgversion) || pkgversion(Functors) < v"0.5"

Functors.@functor Mul
Functors.@functor Add
Functors.@functor MulAdd
Functors.@functor AddMul
Functors.@functor InvMul
Functors.@functor Subtract
Functors.@functor InvMulAdd
Functors.@functor InvAddMul

end # Functors < v"0.5"

end # module AffineMapsFunctorsExt
