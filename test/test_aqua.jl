# This file is a part of AffineMaps.jl, licensed under the MIT License (MIT).

import Test
import Aqua
import AffineMaps

Test.@testset "Aqua tests" begin
    Aqua.test_all(
        AffineMaps,
        ambiguities = true,
        project_toml_formatting = VERSION≥v"1.7"
    )
end # testset
