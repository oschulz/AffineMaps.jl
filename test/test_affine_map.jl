# This file is a part of AffineMaps.jl, licensed under the MIT License (MIT).

using AffineMaps
using Test

using InverseFunctions
using ChangesOfVariables
import ForwardDiff

include("getjacobian.jl")

@testset "AffineMap" begin
    n = 5

    A_scalar = 3.3
    A_scalar_c = Complex(3.3, 1.2)
    A_mat = randn(n, n)
    A_mat_c = Complex.(randn(n, n), randn(n, n))

    b_scalar = 0.7
    b_scalar_c = 0.7
    b_vec = rand(n)
    b_vec_c = Complex.(rand(n), rand(n))

    x_scalar = 0.7
    x_scalar_c = Complex(0.7, 0.3)
    x_vec = rand(n)
    x_vec_c = Complex.(rand(n), rand(n))

    for (A, b, x) in [
        (A_scalar, b_scalar, x_scalar),
        (A_scalar, b_vec, x_vec),
        (A_mat, b_vec, x_vec),

        (A_scalar, b_scalar, x_scalar_c),
        (A_scalar_c, b_scalar, x_scalar_c),
        (A_scalar, b_scalar_c, x_scalar_c),

        (A_scalar, b_vec, x_vec_c),
        (A_scalar_c, b_vec, x_vec_c),
         (A_scalar, b_vec_c, x_vec_c),

        (A_mat, b_vec, x_vec_c),
        (A_mat_c, b_vec, x_vec_c),
        (A_mat, b_vec_c, x_vec_c),
    ]
        for (f, y) in [
            (Mul(A), A * x),
            (Add(b), x + b),
            (MulAdd(A, b), A * x + b),
            (AddMul(b, A), A * (x + b)),
            (InvMul(A), A \ x),
            (Subtract(b), x - b),
            (InvMulAdd(A, b), A \ ( x - b)),
            (InvAddMul(b, A), A \ x - b),
        ]
            global g_state = (;A, b, x, f, y)
            @test @inferred(f(x)) ≈ y
            InverseFunctions.test_inverse(f, x)
            ChangesOfVariables.test_with_logabsdet_jacobian(f, x, getjacobian)
        end
    end

    @test_throws ArgumentError with_logabsdet_jacobian(Mul(A_mat_c), x_vec)
end

nothing
