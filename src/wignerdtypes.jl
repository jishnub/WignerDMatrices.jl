twoj(j) = twice(Integer,j)
function twojp1(j)
	t = twoj(j)
	t + one(t)
end

function filledelements(j)
	N = twojp1(j)
	Np = (N+1) >> 1
	return (1 + N - Np)*Np
end

struct WignerdMatrixContainer{T} <: AbstractMatrix{T}
	j :: HalfInt
	v :: Vector{T}

	function WignerdMatrixContainer{T}(j::HalfInt, v::Vector{T}) where {T}
		
		j >= zero(j) || throw(ArgumentError("j must be ≥ 0"))

		if length(v) < filledelements(j)
			throw(ArgumentError("array provided is not large enough"*
			" to store all values. Need an array of "*
			"length $(filledelements(j)) whereas the array provided "*
			"has length = $(length(v))"))
		end

		new{T}(j, v)
	end
end

function WignerdMatrixContainer{T}(::UndefInitializer, j) where {T}
	jh = HalfInt(j)
	N = filledelements(jh)
	v = Vector{T}(undef, N)
	WignerdMatrixContainer{T}(jh, v)
end

Base.size(d::WignerdMatrixContainer) = (twojp1(d.j), twojp1(d.j))
Base.fill!(d::WignerdMatrixContainer, val) = fill!(d.v, val)
Base.IndexStyle(::Type{<:WignerdMatrixContainer}) = IndexCartesian()

function Base.:(==)(d1::WignerdMatrixContainer, d2::WignerdMatrixContainer)
	d1.j == d2.j || return false
	N = filledelements(d1.j)
	for i = 1:N
		d1.v[i] == d2.v[i] || return false
	end
	return true
end

function flatind(j, m_ind, n_ind)
	# This is m-major, so m increases faster than n
	# Store only the left triangular quadrant
	N = twojp1(j)
	indskip = (2 + N - n_ind)*(n_ind - 1)
	indskip + (m_ind - n_ind) + 1
end

function flatind_phase(j, m_ind, n_ind)
	N = twojp1(j)
	Np = (N + 1) >> 1

	neg = (-1)^(m_ind - n_ind)
	phase = one(neg)

	# The left quadrant (m,n <= 0, n <= m <= -n) is stored
	# We evaluate the others using symmetry relations
	if n_ind > Np && N - n_ind < m_ind <= n_ind 
		# right
		# (m, n) → (-m, -n)
		n_ind, m_ind = N - n_ind + 1, N - m_ind + 1
		phase = neg
	elseif m_ind < n_ind
		# top
		# (m, n) → (n, m)
		m_ind, n_ind = n_ind, m_ind
		phase = neg
	elseif n_ind < m_ind && m_ind + n_ind > N + 1
		# bottom
		# (m, n) → (-n ,-m)
		m_ind, n_ind = N - n_ind + 1, N - m_ind + 1
	end
	
	ind = flatind(j, m_ind, n_ind)
	
	return ind, phase
end

@propagate_inbounds function Base.getindex(d::WignerdMatrixContainer, m_ind::Int, n_ind::Int, I...)
	@boundscheck checkbounds(d, m_ind, n_ind, I...)
	ind, phase = flatind_phase(d.j, m_ind, n_ind)
	d.v[ind]*phase
end

@propagate_inbounds function Base.setindex!(d::WignerdMatrixContainer, val, m_ind::Int, n_ind::Int, I...)
	@boundscheck checkbounds(d, m_ind, n_ind, I...)
	ind, phase = flatind_phase(d.j, m_ind, n_ind)
	d.v[ind] = val*phase
	d
end

abstract type AbstractWignerMatrix{T} <: AbstractHalfIntegerMatrix{T} end

Base.size(w::AbstractWignerMatrix) = size(parent(w))
Base.axes(w::AbstractWignerMatrix) = axes(parent(w))
Base.axes(w::AbstractWignerMatrix, d) = axes(parent(w), d)

"""
	WignerdMatrix(j, beta::Real)

Compute the Wigner d-matrix in the z-y-z convention for the 
angular momentum `j` and the angle `beta`. 
The resulting matrix has a size of `(2j+1, 2j+1)` and axes `(-j:j, -j:j)`. 
The angular momentum `j` needs to be either an `Integer` or a half-integer.
The matrix elements are real in this convention.

# Examples

```jldoctest
julia> d = WignerdMatrix(0.5, π)
2×2 WignerdMatrix{Float64} for j = 1/2 and beta = π with indices -1/2:1/2×-1/2:1/2:
  0.0  1.0
 -1.0  0.0

julia> d[0.5, 0.5]
0.0
```
"""
struct WignerdMatrix{T, B<:Real} <: AbstractWignerMatrix{T}
	beta :: B
	dj :: SpinMatrix{T,WignerdMatrixContainer{T}}
end

Base.parent(d::WignerdMatrix) = d.dj
Base.collect(d::WignerdMatrix) = collect(parent(d))

@propagate_inbounds function Base.getindex(d::WignerdMatrix{T}, m::Real, n::Real) where {T}
	parent(d)[m,n]
end

@propagate_inbounds function Base.getindex(d::WignerdMatrix{T, <:Union{ScaledPi,Irrational{:π}}}, m::Real, n::Real) where {T}
	@boundscheck checkbounds(d, m, n)
	wignerdmatrixelement(T, sphericaldegree(d), m, n, d.beta)
end

#= Linear indexing has to be done in a roundabout way through Cartesian indexing
to avoid undefined references.
=#
@inline function Base.getindex(d::WignerdMatrix, i::Int)
	@boundscheck checkbounds(d, i)
	mn = Tuple(@inbounds eachindex(d)[i])
	@inbounds d[mn...]
end

function Base.:(==)(d1::WignerdMatrix, d2::WignerdMatrix)
	d1.beta == d2.beta && parent(d1.dj) == parent(d2.dj)
end

function Base.:(==)(d1::WignerdMatrix{<:Any,<:ScaledPi}, d2::WignerdMatrix{<:Any,<:ScaledPi})
	(d1.beta == d2.beta && sphericaldegree(d1) == sphericaldegree(d2)) || return false

	j = sphericaldegree(d1)
	for n = -j:zero(j)
		m = nonzerofirstindex(d1.beta, n)
		d1[m,n] == d2[m,n] || return false
	end
	return true
end

"""
	WignerDMatrix(j, α::Real, β::Real, γ::Real)

Construct the Wigner D-matrix in the z-y-z convention for the angular momentum `j` 
and the Euler angles `(α, β, γ)`. The angular momentum `j` needs to be an `Integer` 
or a half-integer. The resulting matrix has a size of `(2j+1, 2j+1)` and axes 
`(-j:j, -j:j)`.

# Examples

```jldoctest
julia> D = WignerDMatrix(0.5, 0, 0, 0)
2×2 WignerDMatrix{Complex{Float64}} for j = 1/2, alpha = 0, beta = 0 and gamma = 0 with indices -1/2:1/2×-1/2:1/2:
 1.0+0.0im  0.0-0.0im
 0.0+0.0im  1.0+0.0im

julia> D[-0.5, 0.5]
-0.0 + 0.0im
```
"""
struct WignerDMatrix{T<:Complex,A<:Real,W<:WignerdMatrix,G<:Real} <: AbstractWignerMatrix{T}
	alpha :: A
	dj :: W
	gamma :: G
end

function WignerDMatrix{T}(α::A, dj::W, γ::G) where {W<:WignerdMatrix,A<:Real,G<:Real,T<:Complex}
	WignerDMatrix{T,A,W,G}(α, dj, γ)
end

function WignerDMatrix(α::A, dj::WignerdMatrix{R}, γ::G) where {R,A<:Real,G<:Real}
	Tphase = Complex{promote_type(HalfInt,promote_type(A,G))}
	T = promote_type(Tphase, R)
	WignerDMatrix{T}(α, dj, γ)
end

Base.parent(D::WignerDMatrix) = D.dj

WignerDphase(m, α::Real) = cis(-m * α)
WignerDphase(m, α::SpecialAngle) = cis_special(-m, α)
WignerDphase(m, α::ZeroRadians) = 1 # return an integer to avoid complex multiplication

@inline function Base.getindex(D::WignerDMatrix{T}, m::Real, n::Real) where {T}
	@boundscheck checkbounds(D, m, n)
	dmn = @inbounds D.dj[m,n]
	val = dmn * WignerDphase(m, D.alpha) * WignerDphase(n, D.gamma)
	convert(T, val)
end

# Linear indexing forces Cartesian indexing
@inline function Base.getindex(D::WignerDMatrix, i::Int)
	@boundscheck checkbounds(D, i)
	mn = Tuple(@inbounds eachindex(D)[i])
	@inbounds D[mn...]
end

function Base.:(==)(d1::WignerDMatrix, d2::WignerDMatrix)
	d1.alpha == d2.alpha && d1.gamma == d2.gamma && d1.dj == d2.dj
end

eulerangles(d::WignerdMatrix) = (zero(d.beta), d.beta, zero(d.beta))
eulerangles(D::WignerDMatrix) = (D.alpha, D.dj.beta, D.gamma)

sphericaldegree(d::WignerdMatrixContainer) = d.j
sphericaldegree(d::SpinMatrix) = d.j
sphericaldegree(d::Union{WignerdMatrix, WignerDMatrix}) = sphericaldegree(d.dj)

Base.IndexStyle(d::Type{<:AbstractWignerMatrix}) = IndexCartesian()

function Base.similar(d::WignerdMatrix, ::Type{T}) where {T}
	jh = HalfInt(sphericaldegree(d))
	w = WignerdMatrixContainer{T}(undef, jh)
	s = SpinMatrix(w, jh)
	WignerdMatrix(d.beta, s)
end

function Base.similar(D::WignerDMatrix, ::Type{T}) where {T}
	d = similar(D.dj)
	WignerDMatrix{T}(D.alpha, d, D.gamma)
end

HalfIntegerArrays.SpinMatrix(d::AbstractWignerMatrix) = SpinMatrix(collect(d))

Base.dataids(d::WignerdMatrixContainer) = Base.dataids(d.v)
Base.dataids(d::AbstractWignerMatrix) = Base.dataids(parent(d))
function Broadcast.broadcast_unalias(dest::AbstractWignerMatrix, src::AbstractWignerMatrix)
	parent(dest) === parent(src) ? src : Broadcast.unalias(dest, src)
end

# Display

function Base.showarg(io::IO, d::WignerdMatrix, toplevel)
	j = sphericaldegree(d)
	if toplevel
		print(io,"WignerdMatrix{$(eltype(d))} for j = ", j, " and beta = ", d.beta)
	else
		print(io,"WignerdMatrix{$(eltype(d))}(",j,", ",d.beta,")")
	end
end
function Base.showarg(io::IO, D::WignerDMatrix, toplevel)
	α, β, γ = eulerangles(D)
	j = sphericaldegree(D)
	if toplevel
		print(io,"WignerDMatrix{$(eltype(D))} for j = ", j,
			", alpha = ", α,", beta = ", β, " and gamma = ", γ)
	else
		print(io,"WignerDMatrix{$(eltype(D))}(", j,
			", ", α, ", ", β, ", ", γ,")")
	end
end

function Base.replace_in_print_matrix(d::WignerdMatrix{<:Any, <:ScaledPi}, m::HalfInteger, n::HalfInteger, s::AbstractString)
    nonzerocheck(d.beta, m, n, sphericaldegree(d)) ? s : Base.replace_with_centered_mark(s)
end

function Base.replace_in_print_matrix(D::WignerDMatrix{<:Any, <:Any, <:WignerdMatrix{<:Any, <:ScaledPi}}, 
	m::HalfInteger, n::HalfInteger, s::AbstractString)

	_,β,_ = eulerangles(D)
	nonzerocheck(β, m, n, sphericaldegree(D)) ? s : Base.replace_with_centered_mark(s)
end