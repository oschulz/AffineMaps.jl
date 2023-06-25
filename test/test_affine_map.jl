# This file is a part of AffineMaps.jl, licensed under the MIT License (MIT).

using AffineMaps
using Test

using LinearAlgebra
using InverseFunctions, ChangesOfVariables
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
    b_mat = randn(n, n)
    b_mat_c = Complex.(randn(n, n), randn(n, n))

    x_scalar = 0.7
    x_scalar_c = Complex(0.7, 0.3)
    x_vec = rand(n)
    x_vec_c = Complex.(rand(n), rand(n))
    x_mat = randn(n, n)
    x_mat_c = Complex.(randn(n, n), randn(n, n))

    for A in [
        A_scalar,
        A_mat,
        A_scalar_c,
        A_mat_c
    ],
    b in [
        b_scalar,
        b_vec,
        b_mat,
        b_scalar_c,
        b_vec_c,
        b_mat_c
    ],
    x in [
        x_scalar,
        x_vec,
        x_mat,
        x_scalar_c,
        x_vec_c,
        x_mat_c
    ]
        for (f, inv_f, y) in [
            (Mul(A), InvMul(A), A * x),
            (Add(b), Subtract(b), x .+ b),
            (MulAdd(A, b), InvMulAdd(A, b), A * x .+ b),
            (AddMul(b, A), InvAddMul(b, A), A * (x .+ b)),
        ]
            if !(A isa AbstractVector && x isa AbstractMatrix)
                @test f isa Function
                @test @inferred(f(x)) ≈ y

                if size(y) == size(x)
                    InverseFunctions.test_inverse(f, x)
                    @test @inferred(inv_f(y)) ≈ x
                    InverseFunctions.test_inverse(inv_f, y)
                    if eltype(A) <: Real && eltype(b) <: Real || eltype(x) <: Complex
                        ChangesOfVariables.test_with_logabsdet_jacobian(f, x, getjacobian)
                        ChangesOfVariables.test_with_logabsdet_jacobian(inv_f, y, getjacobian)
                    end
                end
            end
        end
    end

    @test_throws ArgumentError with_logabsdet_jacobian(Mul(A_mat_c), x_vec)
end

nothing
