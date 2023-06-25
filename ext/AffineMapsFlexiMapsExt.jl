# This file is a part of AffineMaps.jl, licensed under the MIT License (MIT).

module AffineMapsFlexiMapsExt

using FlexiMaps
using AffineMaps

FlexiMaps.islinear(::Mul) = true
FlexiMaps.islinear(::Add) = false
FlexiMaps.islinear(::MulAdd) = false
FlexiMaps.islinear(::AddMul) = false
FlexiMaps.islinear(::InvMul) = true
FlexiMaps.islinear(::Subtract) = false
FlexiMaps.islinear(::InvMulAdd) = false
FlexiMaps.islinear(::InvAddMul) = false

FlexiMaps.isaffine(::AffineMaps.AbstractAffineMap) = true

end # module AffineMapsFlexiMapsExt
