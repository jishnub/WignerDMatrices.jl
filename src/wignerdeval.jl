X(j,n) = sqrt((j+n)*(j-n+1))

function Jy_matrix_zbasis!(j, A::AbstractArray{T}) where {T<:Complex}

	N = twojp1(j)
	length(A) >= N^2 || throw(ArgumentError("array needs to have at least $(N^2) elements"))
	
	Av = reshape(@view(A[1:N^2]),N,N)
	fill!(Av,zero(T))
	h = Hermitian(Av)

    @inbounds for i in 1:N-1
	    Av[i,i+1] = T(0, X(j,-j+i)/2)
	end

	return h
end

function Jy_eigen_nomatch!(j, A)
	Jy = Jy_matrix_zbasis!(j, A)
	_,v = eigen!(Jy, sortby = λ -> round(real(λ)))
	permutedims(v)
end

function eigenvecsJy!(j, A = Matrix{ComplexF64}(undef, twojp1(j), twojp1(j)))
	get!(JyEigenDict, HalfInt(j)) do
		Jy_eigen_nomatch!(j, A)
	end
end

function wignerdmatrixelement(::Type{T}, j::HalfInt, m, n, β::Real, v) where {T<:Real}
	mind = floor(Int, m + j) + 1
	nind = floor(Int, n + j) + 1

	# Compute the sum in the higher accuracy out of T and the individual terms
	Tel = promote_type(eltype(v), Complex{promote_type(HalfInt, typeof(β))})
	Tdcomplex = promote_type(Complex{T},Tel)
	TR_high = typeofreal(Tdcomplex)
	dmn = zero(Tdcomplex)
	
	λ = HalfInt(-j):HalfInt(j)

	@inbounds for μ in eachindex(λ)
		dmn += cis(-TR_high(λ[μ]*β)) * v[μ,mind] * conj(v[μ,nind])
	end

	T(real(dmn))
end

nonzerocheck(::Union{ZeroRadians, TwoNPiRadians}, m, n, j) = m == n
nonzerocheck(::TwoNPlusOnePiRadians, m, n, j) = m == -n
function nonzerocheck(::Piby2Radians, m, n, j)
	!(isodd(floor(Int, j+m)) && iszero(n)) && !(isodd(floor(Int, j+n)) && iszero(m))
end
# Generally this isn't known from the type, assume all are non-zero
nonzerocheck(::Real, m, n, j) = true

nonzerofirstindex(::Union{ZeroRadians, TwoNPiRadians}, n) = n
nonzerofirstindex(::Union{TwoNPlusOnePiRadians}, n) = -n

function wignerdmatrixelement(::Type{T}, j::HalfInt, m, n, β::ZeroRadians, v = nothing) where {T<:Real}
	nonzerocheck(β, m, n, j) ? one(T) : zero(T)
end

function wignerdmatrixelement(::Type{T}, j::HalfInt, m, n, β::Irrational{:π}, v = nothing) where {T<:Real}
	wignerdmatrixelement(T, j, m, n, Pi, v)
end

function wignerdmatrixelement(::Type{T}, j::HalfInt, m, n, x::TwoNPiRadians, v = nothing) where {T<:Real}
	nonzerocheck(x, m, n, j) ? (iseven(twice(j)*x.n) ? one(T) : -one(T)) : zero(T)
end

function wignerdmatrixelement(::Type{T}, j::HalfInt, m, n, x::TwoNPlusOnePiRadians, v = nothing) where {T<:Real}
	nonzerocheck(x, m, n, j) ? (iseven(twice(j)*x.n + floor(Int, j+m)) ? one(T) : -one(T)) : zero(T)
end

function wignerdmatrixelement(::Type{T}, j::HalfInt, m, n, β::Piby2Radians, v) where {T<:Real}
	mind = floor(Int, m + j) + 1
	nind = floor(Int, n + j) + 1

	λ = HalfInt(-j):HalfInt(j)

	# Compute the sum in the higher accuracy out of T and the individual terms
	Tel = promote_type(eltype(v), Complex{promote_type(HalfInt, typeof(β))})
	Tdcomplex = promote_type(Complex{T},Tel)
	TR_high = typeofreal(Tdcomplex)
	dmn = zero(Tdcomplex)
	
	if nonzerocheck(β, m, n, j)
		@inbounds for μ in eachindex(λ)
			dmn += cis_special(TR_high,-λ[μ], β) * v[μ,mind] * conj(v[μ,nind])
		end
	end

	T(real(dmn))
end

typeofreal(::Type{Complex{T}}) where {T} = T

# promote different depending on whether v is specified.
# the result will be real if v is not specified
function wignerdmatrixelement(j::HalfInt, m, n, β::Real, v::AbstractMatrix)
	T = promote_type(eltype(v), Complex{promote_type(HalfInt, typeof(β))})
	wignerdmatrixelement(typeofreal(T), j, m, n, β, v)
end

function wignerdmatrixelement(j::HalfInt, m, n, β::Real, v = nothing)
	T = promote_type(HalfInt, typeof(β))
	wignerdmatrixelement(T, j, m, n, β, v)
end

"""
	wignerdmatrixelement(j::Real, m, n, β::Real, [v = WignerDMatrices.eigenvecsJy!(j)])

Compute the element ``d^j_{m,n}(β)`` of the Wigner d-matrix in the z-y-z convention.
Optionally the matrix `v` that contains the transposed eigenvectors of the angular momentum 
operator `Jy` corresponding to `j` may be provided.

The angular momentum `j` must be either an `Integer` or a half-integer.

# Examples

```jldoctest
julia> dmn = WignerDMatrices.wignerdmatrixelement(1/2, 1/2, 1/2, WignerDMatrices.TwoPi)
-1.0

julia> d = WignerdMatrix(1/2, WignerDMatrices.TwoPi)
2×2 WignerdMatrix{Float64} for j = 1/2 and beta = 2π with indices -1/2:1/2×-1/2:1/2:
 -1.0    ⋅ 
   ⋅   -1.0

julia> d[1/2,1/2] == dmn
true
```
"""
wignerdmatrixelement(j::Real, m, n, β::Real, v = eigenvecsJy!(j)) = wignerdmatrixelement(HalfInt(j), m, n, β, v)

"""
	wignerDmatrixelement(j::Real, m, n, β::Real, [v = WignerDMatrices.eigenvecsJy!(j)])

Compute the element ``D^j_{m,n}(β)`` of the Wigner D-matrix in the z-y-z convention.
Optionally the matrix `v` that contains the transposed eigenvectors of the angular momentum 
operator `Jy` corresponding to `j` may be provided.

The angular momentum `j` must be either an `Integer` or a half-integer.

# Examples

```jldoctest
julia> Dmn = WignerDMatrices.wignerDmatrixelement(1/2, 1/2, 1/2, (0, WignerDMatrices.TwoPi, 0))
-1.0 + 0.0im

julia> D = WignerDMatrix(1/2, 0, WignerDMatrices.TwoPi, 0)
2×2 WignerDMatrix{Complex{Float64}} for j = 1/2, alpha = 0, beta = 2π and gamma = 0 with indices -1/2:1/2×-1/2:1/2:
 -1.0-0.0im       ⋅    
      ⋅      -1.0-0.0im

julia> D[1/2,1/2] == Dmn
true
```
"""
function wignerDmatrixelement(j::Real, m, n, (α, β, γ)::NTuple{3,Real}, v = eigenvecsJy!(j))
	wignerdmatrixelement(j, m, n, β, v) * WignerDphase(m, α) * WignerDphase(n, γ)
end

function djmatrix_fill!(d::SpinMatrix{T,<:WignerdMatrixContainer}, j, β, v) where {T}
	# optimizations
	if iszero(β)
		djmatrix_fill!(d, j, Zero, v)
	elseif isinteger(β/π)
		n = floor(Int, β/π) 
		# type-unstable, but doesn't really matter
		if isodd(n)
			βnew = TwoNPlusOnePiRadians(n >> 1)
		else 
			βnew = TwoNPiRadians(n >> 1)
		end
		djmatrix_fill!(d, j, βnew, v)
	else
		@inbounds for n = -j:zero(j), m = n:-n
			# skip if the element is known to be zero
			nonzerocheck(β, m, n, j) || continue
			d[m,n] = wignerdmatrixelement(T, j, m, n, β, v)
		end
	end

	return d
end

function djmatrix_fill!(d::SpinMatrix{<:Any,<:WignerdMatrixContainer}, j, β::Irrational{:π}, v)
	djmatrix_fill!(d, j, Pi, v)
end

function djmatrix_fill!(d::SpinMatrix{T,<:WignerdMatrixContainer}, j, β::ScaledPi, v) where {T}
	@inbounds for n = -j:zero(j)
		m = nonzerofirstindex(β, n)
		d[m,n] = wignerdmatrixelement(T, j, m, n, β, v)
	end

	return d
end

djmatrix_fill!(d::WignerdMatrix, args...) = djmatrix_fill!(d.dj, args...)
djmatrix_fill!(d::WignerDMatrix, args...) = djmatrix_fill!(d.dj, args...)

function WignerdMatrix!(d::SpinMatrix, β::Real, A)

	fill!(parent(d), zero(eltype(parent(d))))

	j = d.j
	v = eigenvecsJy!(j, A)
	djmatrix_fill!(d, HalfInt(j), β, v)

	return d
end

function WignerdMatrix!(D::WignerDMatrix{T}, β::Real, args...) where {T}
	d = WignerdMatrix!(D.dj, β, args...)
	WignerDMatrix{T}(D.alpha, d, D.gamma)
end

"""
	WignerdMatrix{T}(j, β::Real) where {T<:Real}

Construct the Wigner d-matrix in the z-y-z convention for the angular momentum `j` and the angle `β`, where 
the elements of the matrix will be of the type `T`. 

The matrix elements are expanded in the basis of the eigenfunctions of 
the angular momentum operator `Jy`. 
Specifying the type of the matrix elements does not increase the accuracy 
of the numerical evaluation of eigenvectors, 
as this is set by the LAPACK implementation used.
An increased precision, however, might help in avoiding round-off errors in the sum.

# Examples

```jldoctest
julia> d = WignerdMatrix{BigFloat}(1, 0.4)
3×3 WignerdMatrix{BigFloat} for j = 1 and beta = 0.4 with indices -1:1×-1:1:
  0.96053     0.27536   0.0394695
 -0.27536     0.921061  0.27536
  0.0394695  -0.27536   0.96053

julia> d[1,1]
0.9605304970014426312106516543295635109444272950721034465148443857827477926800066

julia> WignerdMatrix{Float64}(1, 0.4)[1,1] # compare with Float64
0.9605304970014427
```
"""
function WignerdMatrix{T}(j, β::Real) where {T<:Real}
	jh = HalfInt(j)
	w = WignerdMatrixContainer{T}(undef, jh)
	s = SpinMatrix(w, jh)
	d = WignerdMatrix(β, s)

	A = Matrix{ComplexF64}(undef, size(d)...)
	
	WignerdMatrix!(s, β, A)
	return d
end
WignerdMatrix(j, β::T) where {T<:Real} = WignerdMatrix{promote_type(T,Float64)}(j, β)
WignerdMatrix(j, β::SpecialAngle{T}) where {T} = WignerdMatrix{promote_type(T,Float64)}(j, β)

"""
	WignerDMatrix{T}(j, α::Real, β::Real, γ::Real) where {T<:Complex}

Construct the Wigner D-matrix in the z-y-z convention for the angular momentum `j` and the Euler angles 
`(α, β, γ)`, where the matrix elements are of the type `T`.

# Examples

```jldoctest
julia> WignerDMatrix{Complex{BigFloat}}(0.5, 0, 0, 0)
2×2 WignerDMatrix{Complex{BigFloat}} for j = 1/2, alpha = 0, beta = 0 and gamma = 0 with indices -1/2:1/2×-1/2:1/2:
 1.0+0.0im  0.0-0.0im
 0.0+0.0im  1.0+0.0im
```
"""
function WignerDMatrix{Complex{T}}(j, α::Real, β::Real, γ::Real) where {T}
	Td = promote_type(T, Float64)
	d = WignerdMatrix{Td}(j, β)
	WignerDMatrix{Complex{T}}(α, d, γ)
end

function WignerDMatrix(j, α::Real, β::Real, γ::Real)
	d = WignerdMatrix(j, β)
	WignerDMatrix(α, d, γ)
end