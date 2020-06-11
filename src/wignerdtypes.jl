promote_type_phase(T::Type) = promote_type(T,typeof(Integer(zero(T))))

function filledelements(j)
	if isodd(twice(Integer,j))
		h = Integer(j - half(1))
		N = (1 + h)*(2 + h)
	else
		h = Integer(j)
		N = (h + 1)^2
	end
	return N
end

abstract type AbstractWignerMatrix{T} <: AbstractMatrix{T} end

mutable struct WignerdMatrix{T<:Real,B<:Real,V<:AbstractArray{<:Real}} <: AbstractWignerMatrix{T}
	j :: HalfInt
	β :: B
	dj :: V

	function WignerdMatrix{T,B,V}(j::HalfInt, β::B, dj::V) where {T,B<:Real,V<:AbstractArray{<:Real}}
		
		j >= 0 || throw(ArgumentError("j must be ≥ 0"))

		if length(dj) < filledelements(j)
			throw(ArgumentError("array provided is not large enough"*
			" to store all values. Need an array of "*
			"length $(filledelements(j)) whereas the array provided "*
			"has length = $(length(dj))"))
		end

		new{T,B,V}(j, β, dj)
	end
end

function WignerdMatrix(T::Type, j, β::B, dj::V) where {B<:Real,V<:AbstractArray{<:Real}}
	WignerdMatrix{T,B,V}(HalfInt(j), β, dj)
end

function WignerdMatrix(j, β, dj::AbstractArray{R}) where {R<:Real}
	T = promote_type_phase(R)
	WignerdMatrix(T, HalfInt(j), β, dj)
end

function WignerdMatrix{Q}(::UndefInitializer, j, β) where {Q<:Real}
	dj = Vector{Q}(undef, filledelements(j))
	WignerdMatrix(Q, j, β, dj)
end

struct WignerDMatrix{T<:Complex,A<:Real,W<:WignerdMatrix,G<:Real} <: AbstractWignerMatrix{T}
	α :: A
	dj :: W
	γ :: G
end

eulerangles(d::WignerdMatrix) = (zero(d.β), d.β, zero(d.β))
eulerangles(D::WignerDMatrix) = (D.α, D.dj.β, D.γ)

updatebeta!(d::WignerdMatrix,β) = (d.β = β)
updatebeta!(D::WignerDMatrix,β) = updatebeta!(D.dj, β)

function WignerDMatrix(T::Type, α::A, dj::W, γ::G) where {W<:WignerdMatrix,A<:Real,G<:Real}
	WignerDMatrix{T,A,W,G}(α, dj, γ)
end
function WignerDMatrix(α::A, dj::WignerdMatrix{R}, γ::G) where {R<:Real,A<:Real,G<:Real}
	TR = promote_type_phase(R)
	TC = Complex{promote_type(A,G)}
	T = promote_type(TC,TR)
	WignerDMatrix{T,A,WignerdMatrix{R},G}(α, dj, γ)
end

function Base.:(==)(d1::WignerdMatrix, d2::WignerdMatrix)
	d1.j == d2.j && d1.β == d2.β && d1.dj == d2.dj
end
function Base.:(==)(d1::WignerDMatrix, d2::WignerDMatrix)
	d1.α == d2.α && d1.γ == d2.γ && d1.dj == d2.dj
end
function Base.:(==)(d1::WignerdMatrix, d2::WignerDMatrix)
	iszero(d2.α) && iszero(d2.γ) && d1 == d2.dj
end
function Base.:(==)(d1::WignerDMatrix, d2::WignerdMatrix)
	d2 == d1
end
function Base.:(==)(d1::AbstractWignerMatrix, d2::AbstractArray)
	collect(d1) == d2
end
function Base.:(==)(d1::AbstractArray, d2::AbstractWignerMatrix)
	d2 == d1
end

@inline function flatind(j, m, n)
	# This is m-major, so m increases faster than n
	# Store only the left triangular quadrant
	@boundscheck abs(m) > j &&
	throw(ArgumentError("m=$m does not satisfy abs(m) ⩽ j=$j"))
	@boundscheck -j <= n <= 0 ||
	throw(ArgumentError("n=$n does not satisfy -j=$(-j) ⩽ n ⩽ 0"))

	indskip = (2 + Integer(j - n) )*Integer(j + n)
	indskip + Integer(m - n) + 1
end
@inline @Base.propagate_inbounds function flatind(d::WignerdMatrix, m, n)
	flatind(d.j, m, n)
end

@Base.propagate_inbounds function flatind_phase(d::WignerdMatrix, m, n)
	m,n = promote(m,n)
	neg = (-1)^(m-n)
	phase = one(neg)
	# The left is stored
	# We evaluate the other parts using the correct phases
	if n > 0 && abs(m) <= n
		# right
		m,n = -m,-n
		phase = neg
	elseif m < -abs(n)
		# top
		m,n = n,m
		phase = neg
	elseif m > abs(n)
		# bottom
		m,n = -n,-m
	end
	ind = flatind(d, m, n)
	ind,phase
end

function indicescompatible(j,m,n)
	isinteger(j+m) && isinteger(j+n)
end

@inline @Base.propagate_inbounds function Base.getindex(d::WignerDMatrix{T}, m::HalfInteger, n::HalfInteger) where {T}
	val = d.dj[m,n]*cis(-(m*d.α + n*d.γ))
	convert(T, val)
end
@inline @Base.propagate_inbounds function Base.getindex(d::WignerDMatrix{T,<:SpecialAngles}, m::HalfInteger, n::HalfInteger) where {T}
	val = d.dj[m,n]*cis_special(-m,d.α)*cis(-n*d.γ)
	convert(T, val)
end
@inline @Base.propagate_inbounds function Base.getindex(d::WignerDMatrix{T,<:Real,<:WignerdMatrix,G}, 
	m::HalfInteger, n::HalfInteger) where {T,G<:SpecialAngles}

	val = d.dj[m,n]*cis(-m*d.α)*cis_special(-n,d.γ)
	convert(T, val)
end
@inline @Base.propagate_inbounds function Base.getindex(d::WignerDMatrix{T,A,<:WignerdMatrix,G}, 
	m::HalfInteger, n::HalfInteger) where {T,A<:SpecialAngles,G<:SpecialAngles}

	val = d.dj[m,n]*cis_special(-m,d.α)*cis_special(-n,d.γ)
	convert(T, val)
end

@inline Base.@propagate_inbounds function Base.getindex(d::WignerdMatrix{T}, m::HalfInteger, n::HalfInteger) where {T}
	@boundscheck abs(m) <= d.j || throw(ArgumentError("j ⩽ m ⩽ j not satisfied for j = $(d.j) and m = $m"))
	@boundscheck abs(n) <= d.j || throw(ArgumentError("j ⩽ n ⩽ j not satisfied for j = $(d.j) and m = $n"))
	indicescompatible(d.j,m,n) || throw(ArgumentError("d-matrix with j = $(d.j) can not be indexed with m = $m and n = $n"))
	ind,phase = flatind_phase(d,m,n)
	val = d.dj[ind] * phase
	convert(T,val)
end

@inline @Base.propagate_inbounds function Base.getindex(d::AbstractWignerMatrix, m, n)
	getindex(d,HalfInt(m),HalfInt(n))
end

@inline Base.@propagate_inbounds function Base.getindex(d::WignerdMatrix{T}, ind::Integer) where {T}
	val = d.dj[ind]
	convert(T, val)
end

@inline Base.@propagate_inbounds function Base.setindex!(d::WignerDMatrix, val, m::HalfInteger, n::HalfInteger)
	v = val * cis(m*d.α + n*d.γ)
	d.dj[m,n] = real(v)
end
@inline Base.@propagate_inbounds function Base.setindex!(d::WignerdMatrix, val::Real, m::HalfInteger, n::HalfInteger)
	@boundscheck abs(m) <= d.j || throw(ArgumentError("j ⩽ m ⩽ j not satisfied for j = $(d.j) and m = $m"))
	@boundscheck abs(n) <= d.j || throw(ArgumentError("j ⩽ n ⩽ j not satisfied for j = $(d.j) and m = $n"))
	indicescompatible(d.j,m,n) || throw(ArgumentError("d-matrix with j = $(d.j) can not be indexed with m = $m and n = $n"))
	ind, phase = flatind_phase(d,m,n)
	d.dj[ind] = phase * val
end

@inline Base.@propagate_inbounds function Base.setindex!(d::AbstractWignerMatrix, val, m, n)
	setindex!(d,val,HalfInt(m),HalfInt(n))
end

@inline Base.@propagate_inbounds function Base.setindex!(d::WignerdMatrix, val::Real, ind::Integer)
	d.dj[ind] = val
end

sphericaldegree(d::WignerdMatrix) = d.j
sphericaldegree(d::WignerDMatrix) = sphericaldegree(d.dj)

twoj(j) = twice(Integer,j)
function twojp1(j)
	t = twoj(j)
	t + one(t)
end
twojp1(d::AbstractWignerMatrix) = twojp1(sphericaldegree(d))

function Base.axes(d::WignerdMatrix)
	j = d.j
	(-j:j,-j:j)
end
Base.axes(d::WignerDMatrix) = axes(d.dj)
Base.axes(d::AbstractWignerMatrix, dim::Integer) = axes(d)[dim]

Base.size(d::AbstractWignerMatrix) = map(length,axes(d))
Base.size(d::AbstractWignerMatrix, dim::Integer) = length(axes(d)[dim])

function Base.collect(d::AbstractWignerMatrix{T}) where {T}
	dfull = zeros(T, size(d))
	for (indn,n) in enumerate(axes(d,2)), (indm,m) in enumerate(axes(d,1))
		dfull[indm,indn] = d[m,n]
	end
	dfull
end

Base.similar(d::WignerdMatrix{T}) where {T} = WignerdMatrix(T, d.j, d.β, similar(d.dj))
Base.similar(d::WignerdMatrix, ::Type{T}) where {T} = WignerdMatrix(T, d.j, d.β, similar(d.dj))
Base.similar(D::WignerDMatrix{T}) where {T} = WignerDMatrix(T, D.α, similar(D.dj), D.γ)
Base.similar(D::WignerDMatrix, ::Type{T}) where {T} = WignerDMatrix(T, D.α, similar(D.dj), D.γ)

# Display

function Base.summary(io::IO, d::WignerdMatrix)
	print(io,"Wigner d-matrix with j = ",d.j," and β = ",d.β)
end
function Base.summary(io::IO, D::WignerDMatrix)
	α, β, γ = eulerangles(D)
	j = sphericaldegree(D)
	print(io,"Wigner D-matrix with j = ",j,
		", with α = ",α,", β = ",β," and γ = ",γ)
end
function Base.show(io::IO, d::AbstractWignerMatrix)
	summary(io, d)
	println()
	show(io,collect(d))
end
function Base.show(io::IO,m::MIME"text/plain", d::AbstractWignerMatrix)
	summary(io, d)
	println()
	show(io,m,collect(d))
end