# This file is a part of AffineMaps.jl, licensed under the MIT License (MIT).

module AffineMapsAdaptExt

using Adapt
using AffineMaps

Adapt.adapt_structure(target, f::Mul) = Mul(Adapt.adapt(target, f.A))
Adapt.adapt_structure(target, f::InvMul) = InvMul(Adapt.adapt(target, f.A))

Adapt.adapt_structure(target, f::Add) = Add(Adapt.adapt(target, f.b))
Adapt.adapt_structure(target, f::Subtract) = Subtract(Adapt.adapt(target, f.b))

Adapt.adapt_structure(target, f::MulAdd) = MulAdd(Adapt.adapt(target, f.A), Adapt.adapt(target, f.b))
Adapt.adapt_structure(target, f::InvMulAdd) = InvMulAdd(Adapt.adapt(target, f.A), Adapt.adapt(target, f.b))

Adapt.adapt_structure(target, f::AddMul) = AddMul(Adapt.adapt(target, f.b), Adapt.adapt(target, f.A))
Adapt.adapt_structure(target, f::InvAddMul) = InvAddMul(Adapt.adapt(target, f.b), Adapt.adapt(target, f.A))

end # module AffineMapsAdaptExt
