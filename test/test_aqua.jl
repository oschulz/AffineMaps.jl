# This file is a part of AffineTransformations.jl, licensed under the MIT License (MIT).

import Test
import Aqua
import AffineTransformations

Test.@testset "Aqua tests" begin
    Aqua.test_all(AffineTransformations)
end # testset
