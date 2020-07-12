LinearAlgebra.det(d::Union{WignerDMatrix, WignerdMatrix}) = one(eltype(d))

function LinearAlgebra.tr(d::WignerdMatrix{T}) where {T}
	j = sphericaldegree(d)
	N = twojp1(j)
	iszero(d.beta) && return T(N)
	ωby2 = d.beta/2
	if iszero(sin(ωby2))
		# check if odd or even multiple of π
		if isone(cos(ωby2))
			# even multiple of π
			x = N
		else
			# odd multiple of π
			x = (-1)^(N+1)*N
		end
	else
		x = sin(N*ωby2)/sin(ωby2)
	end
	T(x)
end
function LinearAlgebra.tr(d::WignerdMatrix{T,<:ZeroRadians}) where {T}
	j = sphericaldegree(d)
	T(twojp1(j))
end
function LinearAlgebra.tr(d::WignerdMatrix{T, <:TwoNPiRadians}) where {T}
	j = sphericaldegree(d)
	N = twojp1(j)
	t = wignerdmatrixelement(T, j, -j, -j, d.beta)
	N*t
end
function LinearAlgebra.tr(d::WignerdMatrix{T,<:Union{Irrational{:π}, TwoNPlusOnePiRadians}}) where {T}
	j = sphericaldegree(d)
	if isinteger(j)
		res = wignerdmatrixelement(T, j, 0, 0, d.beta)
	else
		res = zero(T)
	end
	return res
end

function LinearAlgebra.tr(D::WignerDMatrix{T}) where {T}
	j = sphericaldegree(D)
	α,β,γ = eulerangles(D)
	cosωby2 = cos(β/2)cos((α+γ)/2)
	if isone(cosωby2)
		return T(twojp1(j))
	elseif isone(-cosωby2)
		N = twojp1(j)
		return (-1)^(N+1) * N
	end
	ωby2 = acos(cosωby2)
	x = sin((2j+1)ωby2)/sin(ωby2)
	T(x)
end
function LinearAlgebra.tr(D::WignerDMatrix{T, <:Any, 
	WignerdMatrix{<:Any, <:Union{Irrational{:π}, TwoNPlusOnePiRadians}}}) where {T}

	j = sphericaldegree(D)
	if isinteger(j)
		res = Complex(tr(parent(D)))
	else
		res = zero(T)
	end
	return res
end

function LinearAlgebra.inv(D::WignerDMatrix)
	j = sphericaldegree(D)
	α, β, γ = eulerangles(D)
	WignerDMatrix(j, -γ, -β, -α)
end

function LinearAlgebra.inv(D::WignerdMatrix)
	j = sphericaldegree(D)
	_, β, _ = eulerangles(D)
	WignerdMatrix(j, -β)
end

function Base.:(*)(d1::WignerdMatrix, d2::WignerdMatrix)
	sphericaldegree(d1) == sphericaldegree(d2) || 
		throw(ArgumentError("d-matrices with j = $(sphericaldegree(d1)) "*
			"and j = $(sphericaldegree(d2)) can not be multiplied"))
	j = sphericaldegree(d1)

	WignerdMatrix(j, d1.beta + d2.beta)
end

function Base.:(*)(D1::WignerDMatrix, D2::WignerDMatrix)
	SpinMatrix(collect(D1) * collect(D2))
end

function Base.:(*)(D1::WignerDMatrix, d2::WignerdMatrix)
	D1 * WignerDMatrix(0, d2, 0)
end
function Base.:(*)(d1::WignerdMatrix, D2::WignerDMatrix)
	WignerDMatrix(0, d1, 0) * D2
end

for f in [:issymmetric, :ishermitian, :isdiag, :isbanded]
	@eval LinearAlgebra.$f(d::AbstractWignerMatrix) = LinearAlgebra.$f(parent(d))
end

function LinearAlgebra.diag(d::AbstractWignerMatrix)
	N = size(d,1)
	v = Vector{eltype(d)}(undef, N)
	j = sphericaldegree(d)
	@inbounds for (i,m) in enumerate(-j:j)
		v[i] = d[m,m]  
	end
	return v
end

function LinearAlgebra.diag(d::WignerdMatrix{T,<:ZeroRadians}) where {T}
	ones(T, size(d,1))
end

function LinearAlgebra.diag(d::WignerdMatrix{T, <:TwoNPiRadians}) where {T}
	j = sphericaldegree(d)
	t = wignerdmatrixelement(T, j, -j, -j, d.beta)
	fill(t, size(d,1))
end

function LinearAlgebra.diag(d::WignerdMatrix{T,
	<:Union{Irrational{:π}, TwoNPlusOnePiRadians}}) where {T}

	N = size(d,1)
	v = zeros(T, N)
	j = sphericaldegree(d)
	if isodd(N) # integer j
		v[div(N,2) + 1] = wignerdmatrixelement(T, j, 0, 0, d.beta)
	end
	return v
end

function LinearAlgebra.diag(d::WignerDMatrix{T, <:Any, 
	<:WignerdMatrix{<:Any, <:Union{Irrational{:π}, TwoNPlusOnePiRadians}}}) where {T}

	N = size(d,1)
	v = zeros(T, N)
	j = sphericaldegree(d)
	_,β,_ = eulerangles(d)
	if isodd(N) # integer j
		v[div(N,2) + 1] = wignerdmatrixelement(typeofreal(T), j, 0, 0, β)
	end
	return v
end