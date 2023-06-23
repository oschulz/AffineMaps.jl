# This file is a part of AffineMaps.jl, licensed under the MIT License (MIT).

import Test

Test.@testset "Package AffineMaps" begin
    include("test_aqua.jl")
    include("test_affine_map.jl")
    include("test_docs.jl")
    isempty(Test.detect_ambiguities(AffineMaps))
end # testset
