# Product
for T in [:AbstractVector, :AbstractMatrix]
	@eval Base.:(*)(d::AbstractWignerMatrix, v::$T) = collect(d)*v
	@eval Base.:(*)(v::$T, d::AbstractWignerMatrix) = v * collect(d)
end
for T in [:Transpose, :Adjoint]
	@eval Base.:(*)(v::$T{U,<:AbstractVector{U}}, d::AbstractWignerMatrix) where {U} = v * collect(d)
end

# Symmetry
for T in [:issymmetric, :ishermitian]
	@eval LinearAlgebra.$T(::AbstractWignerMatrix) = false
end

LinearAlgebra.det(::AbstractWignerMatrix) = one(Float64)

function LinearAlgebra.tr(D::WignerDMatrix{T}) where {T}
	j = sphericaldegree(D)
	α,β,γ = eulerangles(D)
	ωby2 = acos(cos(β/2)cos((α+γ)/2))
	x = sin((2j+1)ωby2)/sin(ωby2)
	iszero(ωby2) ? oftype(x,twojp1(1)) : x
end
function LinearAlgebra.tr(d::WignerdMatrix)
	j = sphericaldegree(d)
	ωby2 = d.β/2
	x = sin((2j+1)ωby2)/sin(ωby2)
	iszero(ωby2) ? oftype(x,twojp1(j)) : x
end
function LinearAlgebra.tr(d::WignerdMatrix{<:Real,ZeroRadians})
	j = sphericaldegree(d)
	float(twojp1(j))
end