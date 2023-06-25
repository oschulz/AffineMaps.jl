# This file is a part of AffineMaps.jl, licensed under the MIT License (MIT).

"""
    AffineMaps

Implements affine transformations in a modular way.
"""
module AffineMaps

using LinearAlgebra

include("affine_map.jl")

@static if !isdefined(Base, :get_extension)
    include("../ext/AffineMapsAdaptExt.jl")
    include("../ext/AffineMapsChangesOfVariablesExt.jl")
    # FlexiMaps supports Julia >= v1.9 only
    include("../ext/AffineMapsFunctorsExt.jl")
    include("../ext/AffineMapsInverseFunctionsExt.jl")
end

end # module
