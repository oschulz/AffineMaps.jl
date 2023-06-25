# This file is a part of AffineMaps.jl, licensed under the MIT License (MIT).

module AffineMapsInverseFunctionsExt

using InverseFunctions
using AffineMaps


InverseFunctions.inverse(f::Mul) = InvMul(f.A)
InverseFunctions.inverse(f::InvMul) = Mul(f.A)

InverseFunctions.inverse(f::Add) = Subtract(f.b)
InverseFunctions.inverse(f::Subtract) = Add(f.b)

InverseFunctions.inverse(f::MulAdd) = InvMulAdd(f.A, f.b)
InverseFunctions.inverse(f::InvMulAdd) = MulAdd(f.A, f.b)

InverseFunctions.inverse(f::AddMul) = InvAddMul(f.b, f.A)
InverseFunctions.inverse(f::InvAddMul) = AddMul(f.b, f.A)


end # module AffineMapsInverseFunctionsExt
