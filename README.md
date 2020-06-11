# WignerDMatrices

[![Build Status](https://travis-ci.com/jishnub/WignerDMatrices.jl.svg?branch=master)](https://travis-ci.com/jishnub/WignerDMatrices.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/jishnub/WignerDMatrices.jl?svg=true)](https://ci.appveyor.com/project/jishnub/WignerDMatrices-jl)
[![Codecov](https://codecov.io/gh/jishnub/WignerDMatrices.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/jishnub/WignerDMatrices.jl)
[![Coveralls](https://coveralls.io/repos/github/jishnub/WignerDMatrices.jl/badge.svg?branch=master)](https://coveralls.io/github/jishnub/WignerDMatrices.jl?branch=master)

Wigner D and d matrices following [Feng et al. (2015)](https://arxiv.org/abs/1507.04535) and the phase convention used by Varshalovich [1988]. The D matrix is defined as 

`Djmn(α, β, γ) = djmn(β) * exp(-i*(m*α + n*γ))`

This carries out a passive rotation by the Euler angles `(α, β, γ)` in the `z-y-z` convention. The d-matrix is purely real in this convention. The d-matrix is evaluated by explicitly diagonalizing the angular momentum operator `Jy` in the basis of the normalized eigenvectors of `Jz` for a specified degree `j`. The result is expected to be limited by machine accuracy.

# Creating the matrices

Compute the Wigner d-matrix using the signature `WignerdMatrix(j, β)` and the D-matrix using `WignerdMatrix(j, α, β, γ)`. Here the angular momentum `j` may either be an integer or a half-integer.

```julia
julia> d = WignerdMatrix(2, π/3)
Wigner d-matrix with j = 2 and β = 1.0471975511965976
5×5 Array{Float64,2}:
  0.5625     0.649519      0.459279   0.216506     0.0625
 -0.649519   5.82867e-16   0.53033    0.5          0.216506
  0.459279  -0.53033      -0.125      0.53033      0.459279
 -0.216506   0.5          -0.53033    5.82867e-16  0.649519
  0.0625    -0.216506      0.459279  -0.649519     0.5625

julia> D = WignerDMatrix(1/2, 0, π/3, 0.5)
Wigner D-matrix with j = 1/2, with α = 0, β = 1.0471975511965976 and γ = 0.5
2×2 Array{Complex{Float64},2}:
  0.839103+0.214258im  0.484456-0.123702im
 -0.484456-0.123702im  0.839103-0.214258im
```

Optionally the element type of the matrix may be specified as the first argument.

The value of `β` may be updated using the mutating function `WignerdMatrix!`.

```julia
julia> WignerdMatrix!(d, π/10)
Wigner d-matrix with j = 2 and β = 0.3141592653589793
5×5 Array{Float64,2}:
  0.951655      0.301455     0.0584764   0.00756218  0.000598866
 -0.301455      0.880037     0.359943    0.0710198   0.00756218
  0.0584764    -0.359943     0.856763    0.359943    0.0584764
 -0.00756218    0.0710198   -0.359943    0.880037    0.301455
  0.000598866  -0.00756218   0.0584764  -0.301455    0.951655
```

While the values of `α` and `γ` are stored in a `WignerDMatrix`, that of `β` is stored in the underlying `WignerdMatrix`. The way to obtain all the angles is 

```julia
julia> D = WignerDMatrix(1/2, 0, π/3, 0.5);

julia> WignerDMatrices.eulerangles(D)
(0, 1.0471975511965976, 0.5)
```

# Indexing

The matrices may be indexed using appropriate `(m,n)` pairs. These must be integers if `j` is an integer, or half-integers if `j` is a half-integer as well.

```julia
julia> d[1, 1]
0.8800367553350505

julia> D[-1/2, 1/2]
0.4844562108553222 - 0.12370197962726143im
```

The structs are defined to respect the symmetries of the d-matrix.

# Special Angles

The package provides the three special angles `ZeroRadians`, `Piby2Radians` and `PiRadians`. Using these might lead to certain elements of the d-matrix being explicitly evaluated to zero.

```julia
julia> WignerdMatrix(1, WignerDMatrices.Piby2Radians())
Wigner d-matrix with j = 1 and β = Piby2Radians()
3×3 Array{Float64,2}:
  0.5        0.707107  0.5
 -0.707107   0.0       0.707107
  0.5       -0.707107  0.5

julia> WignerdMatrix(1, π/2)
Wigner d-matrix with j = 1 and β = 1.5707963267948966
3×3 Array{Float64,2}:
  0.5        0.707107     0.5
 -0.707107   6.12323e-17  0.707107
  0.5       -0.707107     0.5

julia> WignerdMatrix(1, WignerDMatrices.PiRadians())
Wigner d-matrix with j = 1 and β = PiRadians()
3×3 Array{Float64,2}:
 0.0  -0.0   1.0
 0.0  -1.0  -0.0
 1.0   0.0   0.0

julia> WignerdMatrix(1, π)
Wigner d-matrix with j = 1 and β = π
3×3 Array{Float64,2}:
 -2.27596e-15   8.65956e-17   1.0
 -8.65956e-17  -1.0           8.65956e-17
  1.0          -8.65956e-17  -2.27596e-15
```

Note that using the mutating function `WignerMatrix!` on a d-matrix defined using a special angle where the new angle is a floating-point number will not work.

# Limitations

 * Broadcasting does not work currently.
 * Much of linear-algebra is not implemented for half-integer axes