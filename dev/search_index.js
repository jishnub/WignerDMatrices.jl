var documenterSearchIndex = {"docs":
[{"location":"","page":"Reference","title":"Reference","text":"CurrentModule = WignerDMatrices","category":"page"},{"location":"#WignerDMatrices.jl","page":"Reference","title":"WignerDMatrices.jl","text":"","category":"section"},{"location":"","page":"Reference","title":"Reference","text":"Modules = [WignerDMatrices]","category":"page"},{"location":"#WignerDMatrices.WignerDMatrix","page":"Reference","title":"WignerDMatrices.WignerDMatrix","text":"WignerDMatrix(j, α::Real, β::Real, γ::Real)\n\nConstruct the Wigner D-matrix in the z-y-z convention for the angular momentum j  and the Euler angles (α, β, γ). The angular momentum j needs to be an Integer  or a half-integer. The resulting matrix has a size of (2j+1, 2j+1) and axes  (-j:j, -j:j).\n\nExamples\n\njulia> D = WignerDMatrix(0.5, 0, 0, 0)\n2×2 WignerDMatrix{Complex{Float64}} for j = 1/2, alpha = 0, beta = 0 and gamma = 0 with indices -1/2:1/2×-1/2:1/2:\n 1.0+0.0im  0.0-0.0im\n 0.0+0.0im  1.0+0.0im\n\njulia> D[-0.5, 0.5]\n-0.0 + 0.0im\n\n\n\n\n\n","category":"type"},{"location":"#WignerDMatrices.WignerDMatrix-Union{Tuple{T}, Tuple{Any,Real,Real,Real}} where T","page":"Reference","title":"WignerDMatrices.WignerDMatrix","text":"WignerDMatrix{T}(j, α::Real, β::Real, γ::Real) where {T<:Complex}\n\nConstruct the Wigner D-matrix in the z-y-z convention for the angular momentum j and the Euler angles  (α, β, γ), where the matrix elements are of the type T.\n\nExamples\n\njulia> WignerDMatrix{Complex{BigFloat}}(0.5, 0, 0, 0)\n2×2 WignerDMatrix{Complex{BigFloat}} for j = 1/2, alpha = 0, beta = 0 and gamma = 0 with indices -1/2:1/2×-1/2:1/2:\n 1.0+0.0im  0.0-0.0im\n 0.0+0.0im  1.0+0.0im\n\n\n\n\n\n","category":"method"},{"location":"#WignerDMatrices.WignerdMatrix","page":"Reference","title":"WignerDMatrices.WignerdMatrix","text":"WignerdMatrix(j, beta::Real)\n\nCompute the Wigner d-matrix in the z-y-z convention for the  angular momentum j and the angle beta.  The resulting matrix has a size of (2j+1, 2j+1) and axes (-j:j, -j:j).  The angular momentum j needs to be either an Integer or a half-integer. The matrix elements are real in this convention.\n\nExamples\n\njulia> d = WignerdMatrix(0.5, π)\n2×2 WignerdMatrix{Float64} for j = 1/2 and beta = π with indices -1/2:1/2×-1/2:1/2:\n  0.0  1.0\n -1.0  0.0\n\njulia> d[0.5, 0.5]\n0.0\n\n\n\n\n\n","category":"type"},{"location":"#WignerDMatrices.WignerdMatrix-Union{Tuple{T}, Tuple{Any,Real}} where T<:Real","page":"Reference","title":"WignerDMatrices.WignerdMatrix","text":"WignerdMatrix{T}(j, β::Real) where {T<:Real}\n\nConstruct the Wigner d-matrix in the z-y-z convention for the angular momentum j and the angle β, where  the elements of the matrix will be of the type T. \n\nThe matrix elements are expanded in the basis of the eigenfunctions of  the angular momentum operator Jy.  Specifying the type of the matrix elements does not increase the accuracy  of the numerical evaluation of eigenvectors,  as this is set by the LAPACK implementation used. An increased precision, however, might help in avoiding round-off errors in the sum.\n\nExamples\n\njulia> d = WignerdMatrix{BigFloat}(1, 0.4)\n3×3 WignerdMatrix{BigFloat} for j = 1 and beta = 0.4 with indices -1:1×-1:1:\n  0.96053     0.27536   0.0394695\n -0.27536     0.921061  0.27536\n  0.0394695  -0.27536   0.96053\n\njulia> d[1,1]\n0.9605304970014426312106516543295635109444272950721034465148443857827477926800066\n\njulia> WignerdMatrix{Float64}(1, 0.4)[1,1] # compare with Float64\n0.9605304970014427\n\n\n\n\n\n","category":"method"},{"location":"#WignerDMatrices.wignerDmatrixelement","page":"Reference","title":"WignerDMatrices.wignerDmatrixelement","text":"wignerDmatrixelement(j::Real, m, n, β::Real, [v = WignerDMatrices.eigenvecsJy!(j)])\n\nCompute the element D^j_mn(β) of the Wigner D-matrix in the z-y-z convention. Optionally the matrix v that contains the transposed eigenvectors of the angular momentum  operator Jy corresponding to j may be provided.\n\nThe angular momentum j must be either an Integer or a half-integer.\n\nExamples\n\njulia> Dmn = WignerDMatrices.wignerDmatrixelement(1/2, 1/2, 1/2, (0, WignerDMatrices.TwoPi, 0))\n-1.0 + 0.0im\n\njulia> D = WignerDMatrix(1/2, 0, WignerDMatrices.TwoPi, 0)\n2×2 WignerDMatrix{Complex{Float64}} for j = 1/2, alpha = 0, beta = 2π and gamma = 0 with indices -1/2:1/2×-1/2:1/2:\n -1.0-0.0im       ⋅    \n      ⋅      -1.0-0.0im\n\njulia> D[1/2,1/2] == Dmn\ntrue\n\n\n\n\n\n","category":"function"},{"location":"#WignerDMatrices.wignerdmatrixelement","page":"Reference","title":"WignerDMatrices.wignerdmatrixelement","text":"wignerdmatrixelement(j::Real, m, n, β::Real, [v = WignerDMatrices.eigenvecsJy!(j)])\n\nCompute the element d^j_mn(β) of the Wigner d-matrix in the z-y-z convention. Optionally the matrix v that contains the transposed eigenvectors of the angular momentum  operator Jy corresponding to j may be provided.\n\nThe angular momentum j must be either an Integer or a half-integer.\n\nExamples\n\njulia> dmn = WignerDMatrices.wignerdmatrixelement(1/2, 1/2, 1/2, WignerDMatrices.TwoPi)\n-1.0\n\njulia> d = WignerdMatrix(1/2, WignerDMatrices.TwoPi)\n2×2 WignerdMatrix{Float64} for j = 1/2 and beta = 2π with indices -1/2:1/2×-1/2:1/2:\n -1.0    ⋅ \n   ⋅   -1.0\n\njulia> d[1/2,1/2] == dmn\ntrue\n\n\n\n\n\n","category":"function"}]
}
