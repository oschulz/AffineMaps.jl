# This file is a part of AffineMaps.jl, licensed under the MIT License (MIT).

using AffineMaps
using Test

using LinearAlgebra
using InverseFunctions, ChangesOfVariables
import Adapt, Functors
import ForwardDiff
import FlexiMaps

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

        B = A .+ 1; c = b .+ 1

        @test Mul(A) != Mul(B)
        @test InvMul(A) != InvMul(B)
        @test Add(A) != Add(B)
        @test Subtract(A) != Subtract(B)
        @test MulAdd(A, b) != MulAdd(B, b)
        @test MulAdd(A, b) != MulAdd(A, c)
        @test InvMulAdd(A, b) != InvMulAdd(B, b)
        @test InvMulAdd(A, b) != InvMulAdd(A, c)
        @test AddMul(b, A) != AddMul(c, A)
        @test AddMul(b, A) != AddMul(b, B)
        @test InvAddMul(b, A) != InvAddMul(c, A)
        @test InvAddMul(b, A) != InvAddMul(b, B)

        @test !isapprox(Mul(A), Mul(B))
        @test !isapprox(InvMul(A), InvMul(B))
        @test !isapprox(Add(A), Add(B))
        @test !isapprox(Subtract(A), Subtract(B))
        @test !isapprox(MulAdd(A, b), MulAdd(B, c))
        @test !isapprox(InvMulAdd(A, b), InvMulAdd(B, c))
        @test !isapprox(AddMul(b, A), AddMul(c, B))
        @test !isapprox(InvAddMul(b, A), InvAddMul(c, B))
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
            Ax = A * x
            for (f, inv_f, y, valid) in [
                (Mul(A), InvMul(A), Ax, true),
                (Add(b), Subtract(b), x .+ b, axes(x .+ b) == axes(x)),
                (MulAdd(A, b), InvMulAdd(A, b), Ax .+ b, axes(Ax .+ b) == axes(Ax)),
                (AddMul(b, A), InvAddMul(b, A), A * (x .+ b), axes(x .+ b) == axes(x)),
            ]
                @testset "$(typeof(f)) with x::$(typeof(x))" begin
                    @test f isa Function
                    if !valid
                        @test_throws DimensionMismatch f(x)
                    else
                        @test @inferred(f(x)) ≈ y

                        if size(y) == size(x)
                            InverseFunctions.test_inverse(f, x)
                            @test @inferred(inv_f(y)) ≈ x
                            InverseFunctions.test_inverse(inv_f, y)
                            if (eltype(A) <: Real && eltype(b) <: Real || eltype(x) <: Complex) && !(x isa Matrix)
                                ChangesOfVariables.test_with_logabsdet_jacobian(f, x, getjacobian)
                                ChangesOfVariables.test_with_logabsdet_jacobian(inv_f, y, getjacobian)
                            elseif (eltype(A) <: Real && eltype(b) <: Real || eltype(x) <: Complex) && (x isa AbstractMatrix)
                                for (g, u, v_ref) in ((f, x, y), (inv_f, y, x))
                                    v, ladj = ChangesOfVariables.with_logabsdet_jacobian(g, u)
                                    @test v ≈ v_ref
                                    @test ladj isa AbstractMatrix && size(ladj) == (1, size(u, 2))
                                    @test sum(ladj) ≈ logabsdet(getjacobian(g, u))[1]
                                    if !(b isa AbstractMatrix)
                                        @test all(ladj[1, j] ≈ ChangesOfVariables.with_logabsdet_jacobian(g, u[:, j])[2] for j in axes(u, 2))
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    @testset "UniformScaling operator" begin
        b = rand(n)
        x = rand(n)
        X = rand(n, n)

        for λ in (3.3, 3.3 + 1.2im)
            xl = λ isa Complex ? complex.(x) : x
            for (f, y) in [
                (Mul(λ * I), λ * xl),
                (MulAdd(λ * I, b), λ * xl .+ b),
                (AddMul(b, λ * I), λ * (xl .+ b)),
            ]
                @test @inferred(f(xl)) ≈ y
                InverseFunctions.test_inverse(f, xl)
                ChangesOfVariables.test_with_logabsdet_jacobian(f, xl, getjacobian)
                ChangesOfVariables.test_with_logabsdet_jacobian(inverse(f), y, getjacobian)
            end
        end

        y, ladj = ChangesOfVariables.with_logabsdet_jacobian(Mul(3.3 * I), X)
        @test y ≈ 3.3 * X
        @test ladj ≈ fill(n * log(3.3), 1, n)
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

            params, f_ctor = @inferred Functors.functor(f)
            @test @inferred(f_ctor(params)) == f
            @test Functors.fmap(Array{Float32}, f) isa AffineMaps.AbstractAffineMap
            @test Functors.fmap(Array{Float32}, f) ≈ f
            @test @inferred(Functors.fmap(Array{Float32}, f)(Float32.(x))) isa Vector{Float32}

            @test @inferred(FlexiMaps.isaffine(f)) == true
            @test @inferred(FlexiMaps.islinear(f)) == (f(x + x2) ≈ f(x) + f(x2))
        end
    end
end

nothing
