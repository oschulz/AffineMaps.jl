# This file is a part of AffineMaps.jl, licensed under the MIT License (MIT).

"""
    AffineMaps

Implements affine transformations in a modular way.
"""
module AffineMaps

using LinearAlgebra

using InverseFunctions
using ChangesOfVariables
using FunctionChains

include("affine_step.jl")
include("affine_map.jl")

end # module
