# This file is a part of AffineMaps.jl, licensed under the MIT License (MIT).

import Test
import Aqua
import AffineMaps

Test.@testset "Aqua tests" begin
    Aqua.test_all(
        AffineMaps,
        ambiguities = true
    )
end # testset
