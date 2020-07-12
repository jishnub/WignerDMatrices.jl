# WignerDMatrices

[![Build Status](https://travis-ci.com/jishnub/WignerDMatrices.jl.svg?branch=master)](https://travis-ci.com/jishnub/WignerDMatrices.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/jishnub/WignerDMatrices.jl?svg=true)](https://ci.appveyor.com/project/jishnub/WignerDMatrices-jl)
[![Codecov](https://codecov.io/gh/jishnub/WignerDMatrices.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/jishnub/WignerDMatrices.jl)
[![Coveralls](https://coveralls.io/repos/github/jishnub/WignerDMatrices.jl/badge.svg?branch=master)](https://coveralls.io/github/jishnub/WignerDMatrices.jl?branch=master)

Wigner D and d matrices following [Feng et al. (2015)](https://arxiv.org/abs/1507.04535) and the phase convention used by Varshalovich [1988]. The D matrix is defined as 

`Djmn(α, β, γ) = djmn(β) * exp(-i*(m*α + n*γ))`

This carries out a passive rotation by the Euler angles `(α, β, γ)` in the `z-y-z` convention. The d-matrix is purely real in this convention. The d-matrix is evaluated by explicitly diagonalizing the angular momentum operator `Jy` in the basis of the normalized eigenvectors of `Jz` for a specified degree `j`. The result is expected to be limited by machine accuracy.

# Creating the matrices

Compute the Wigner d-matrix using the signature `WignerdMatrix(j, β)` and the D-matrix using `WignerdMatrix(j, α, β, γ)`, where the angular momentum `j` may either be an integer or a half-integer.

```julia
julia> d = WignerdMatrix(2, π/3)
5×5 WignerdMatrix{Float64} for j = 2 and beta = 1.0471975511965976 with indices -2:2×-2:2:
  0.5625     0.649519      0.459279   0.216506     0.0625
 -0.649519   5.82867e-16   0.53033    0.5          0.216506
  0.459279  -0.53033      -0.125      0.53033      0.459279
 -0.216506   0.5          -0.53033    5.82867e-16  0.649519
  0.0625    -0.216506      0.459279  -0.649519     0.5625

julia> D = WignerDMatrix(1/2, 0, π/3, 0.5)
2×2 WignerDMatrix{Complex{Float64}} for j = 1/2, alpha = 0, beta = 1.0471975511965976 and gamma = 0.5 with indices -1/2:1/2×-1/2:1/2:
  0.839103+0.214258im  0.484456-0.123702im
 -0.484456-0.123702im  0.839103-0.214258im
```

The way to obtain the Euler angles corresponding to the D or the d matrix is using the function `eulerangles`.

```julia
julia> D = WignerDMatrix(1/2, 0, π/3, 0.5);

julia> WignerDMatrices.eulerangles(D)
(0, 1.0471975511965976, 0.5)
```

# Indexing

The matrices may be indexed using appropriate `(m,n)` pairs. These must be integers if `j` is an integer, or half-integers if `j` is a half-integer.

```julia
julia> d = WignerdMatrix(2, π/3);

julia> d[-2,-2]
0.5625000000000003

julia> D = WignerDMatrix(1/2, 0, π/3, 0.5);

julia> D[-1/2, 1/2]
0.4844562108553222 - 0.12370197962726143im
```

The structs are defined to respect the symmetries of the d-matrix. In particular, this means that `d[-m,-n] = (-1)^(m-n) d[m,n]`, `d[n,m] = (-1)^(m-n) d[m,n]` and `d[-n,-m] = d[m,n]`.

# Special Angles

The package provides the special angles `ZeroRadians`, `Piby2Radians`, `PiRadians`, `TwoNPiRadians` and `TwoNPlusOnePiRadians`, and correspondingly the constants `Zero`, `Pi`, `TwoPi`, `ThreePi`, `FourPi` and `Piby2`. Using these might lead to certain elements of the d-matrix being explicitly evaluated to zero, which would lead to faster construction and indexing. The generic types are defined as `TwoNPiRadians(m) = 2mπ` and `TwoNPlusOnePiRadians(m) = (2m + 1)π`. Neither the types nor the constants are exported to avoid conflicts with other packages.

```julia
julia> WignerdMatrix(1, π/2)[0,0]
6.12323399573689e-17

julia> WignerdMatrix(1, WignerDMatrices.Piby2)[0,0]
0.0

julia> WignerdMatrix(1, π)
3×3 WignerdMatrix{Float64} for j = 1 and beta = π with indices -1:1×-1:1:
 0.0   0.0  1.0
 0.0  -1.0  0.0
 1.0   0.0  0.0

julia> WignerdMatrix(1, WignerDMatrices.Pi)
3×3 WignerdMatrix{Float64} for j = 1 and beta = Pi with indices -1:1×-1:1:
  ⋅     ⋅   1.0
  ⋅   -1.0   ⋅ 
 1.0    ⋅    ⋅ 

# Indexing is faster

julia> d1 = WignerdMatrix(40, WignerDMatrices.Zero);

julia> d2 = WignerdMatrix(40, 0);

julia> @btime $d1[1,1];
  4.797 ns (0 allocations: 0 bytes)

julia> @btime $d2[1,1];
  13.249 ns (0 allocations: 0 bytes)
```

# Limitations

 * WignerDMatrix types are not closed on matrix multiplication
 * Much of linear-algebra is not implemented for half-integer axes