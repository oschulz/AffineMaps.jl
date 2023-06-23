# AffineMaps.jl

AffineMaps.jl implements Affine maps. It provides the function object types [`Mul`](@ref), [`Add`](@ref), [`MulAdd`](@ref) and [`AddMul`](@ref), as well as their inverses.

All function objects defined here support the implement/support the APIs of

* [InverseFunctions.jl](https://github.com/JuliaMath/InverseFunctions.jl)
* [ChangesOfVariables.jl](https://github.com/JuliaMath/ChangesOfVariables.jl)
* [Functors.jl](https://github.com/FluxML/Functors.jl)

Example:

```jldoctest example; output = false
using AffineMaps
using LinearAlgebra, InverseFunctions, ChangesOfVariables
A = rand(5, 5)
b = rand(5)
x = rand(5)

f = MulAdd(A, b)
y = f(x)
y ≈ A * x + b

# output

true
```

```jldoctest example; output = false
inverse(f)(y) ≈ x

# output

true
```

```jldoctest example; output = false
y, ladj = with_logabsdet_jacobian(f, x)
y ≈ A * x + b && ladj ≈ logabsdet(A)[1]

# output

true
```
