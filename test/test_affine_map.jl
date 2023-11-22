# This file is a part of AffineMaps.jl, licensed under the MIT License (MIT).

using AffineMaps
using Test

using LinearAlgebra
using InverseFunctions, ChangesOfVariables
import Adapt, Functors
import ForwardDiff

import Pkg
if ("FlexiMaps" in keys(Pkg.project().dependencies))
    # FlexiMaps supports Julia >= v1.9 only.
    import FlexiMaps
end

include("getjacobian.jl")

@testset "AffineMap" begin
    n = 5

    @testset "equality" begin
        A = randn(n, n); A2 = copy(A)
        b = randn(n); b2 = copy(b)

        @test @inferred(Mul(A) == Mul(A2))
        @test @inferred(InvMul(A) == InvMul(A2))
        @test @inferred(Add(A) == Add(A2))
        @test @inferred(Subtract(A) == Subtract(A2))
        @test @inferred(MulAdd(A, b) == MulAdd(A2, b2))
        @test @inferred(InvMulAdd(A, b) == InvMulAdd(A2, b2))
        @test @inferred(AddMul(b, A) == AddMul(b2, A2))
        @test @inferred(InvAddMul(b, A) == InvAddMul(b2, A2))

        @test @inferred(isapprox(Mul(A), Mul(A2); atol = 1e-5))
        @test @inferred(isapprox(InvMul(A), InvMul(A2); atol = 1e-5))
        @test @inferred(isapprox(Add(A), Add(A2); atol = 1e-5))
        @test @inferred(isapprox(Subtract(A), Subtract(A2); atol = 1e-5))
        @test @inferred(isapprox(MulAdd(A, b), MulAdd(A2, b2); atol = 1e-5))
        @test @inferred(isapprox(InvMulAdd(A, b), InvMulAdd(A2, b2); atol = 1e-5))
        @test @inferred(isapprox(AddMul(b, A), AddMul(b2, A2); atol = 1e-5))
        @test @inferred(isapprox(InvAddMul(b, A), InvAddMul(b2, A2); atol = 1e-5))
    end    

    @testset "functionality" begin
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

        @test_throws ArgumentError with_logabsdet_jacobian(Mul(A_mat_c), x_vec)

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
                            #ChangesOfVariables.test_with_logabsdet_jacobian(f, x, getjacobian)
                            #ChangesOfVariables.test_with_logabsdet_jacobian(inv_f, y, getjacobian)
                            @test isapprox(ChangesOfVariables.with_logabsdet_jacobian(f, x)[1], y) && all(isapprox.(ChangesOfVariables.with_logabsdet_jacobian(f, x)[2], logabsdet(getjacobian(f, x))[1]))
                            @test isapprox(ChangesOfVariables.with_logabsdet_jacobian(inv_f, y)[1], x) && all(isapprox.(ChangesOfVariables.with_logabsdet_jacobian(inv_f, y)[2], logabsdet(getjacobian(inv_f, y))[1]))
                            
                        end
                    end
                end
            end
        end
    end

    @testset "Extensions" begin
        A = randn(n, n)
        b = randn(n)
        x = randn(n)
        x2 = randn(n)

        for f in [
                Mul(A),
                InvMul(A),
                Add(b),
                Subtract(b),
                MulAdd(A, b),
                InvMulAdd(A, b),
                AddMul(b, A),
                InvAddMul(b, A)
        ]

            @test @inferred(Adapt.adapt(Array{Float32}, f)) isa AffineMaps.AbstractAffineMap
            @test Adapt.adapt(Array{Float32}, f) ≈ f
            @test @inferred(Adapt.adapt(Array{Float32}, f)(Float32.(x))) isa Vector{Float32}

            @test Functors.fmap(Array{Float32}, f) isa AffineMaps.AbstractAffineMap
            @test Functors.fmap(Array{Float32}, f) ≈ f
            @test @inferred(Functors.fmap(Array{Float32}, f)(Float32.(x))) isa Vector{Float32}

            @static if isdefined(Main, :FlexiMaps)
                @test @inferred(FlexiMaps.isaffine(f)) == true
                @test @inferred(FlexiMaps.islinear(f)) == (f(x + x2) ≈ f(x) + f(x2))
            end
        end
    end
end

nothing
