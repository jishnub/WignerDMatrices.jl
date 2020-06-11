X(j,n) = sqrt((j+n)*(j-n+1))

function Jy_matrix_zbasis!(j, A::AbstractArray{T}) where {T<:Complex}

	N = twojp1(j)
	length(A) >= N^2 || throw(ArgumentError("array isn't long enough to store Jy"))
	
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

function Jy_eigen!(j, A)
	get!(JyEigenDict, HalfInt(j)) do
		Jy_eigen_nomatch!(j, A)
	end
end

function djmatrix_terms(β::Real,λ,v,m,n,j)
	dj_m_n = zero(ComplexF64)

	mind = searchsortedfirst(-j:j, m)
	nind = searchsortedfirst(-j:j, n)

	for μ in axes(λ,1)
		dj_m_n += cis(-λ[μ]*β) * v[μ,mind] * conj(v[μ,nind])
	end

	dj_m_n
end

function djmatrix_terms(β::ZeroRadians,λ,v,m,n,j)
	(m == n) ? one(ComplexF64) : zero(ComplexF64)
end

function djmatrix_terms(β::PiRadians,λ,v,m,n,j)
	(m == -n) ? ComplexF64((-1)^(j+m)) : zero(ComplexF64)
end

function djmatrix_terms(β::Piby2Radians,λ,v,m,n,j)
	dj_m_n = zero(ComplexF64)

	mind = searchsortedfirst(-j:j, m)
	nind = searchsortedfirst(-j:j, n)

	if !(isodd(Integer(j+m)) && n == 0) && !(isodd(Integer(j+n)) && m == 0)
		for μ in axes(λ,1)
			dj_m_n += cis_special(-λ[μ],β) * v[μ,mind] * conj(v[μ,nind])
		end
	end

	dj_m_n
end

function djmatrix_fill!(d::WignerdMatrix,j,β,λ,v)
	for n = -j:zero(j), m = n:-n
		d[m,n] = real(djmatrix_terms(β,λ,v,m,n,j))
	end

	return d
end
djmatrix_fill!(d::WignerDMatrix,args...) = djmatrix_fill!(d.dj,args...)

function WignerdMatrix!(d::AbstractWignerMatrix, β::Real,
	A = Matrix{ComplexF64}(undef,twojp1(d),twojp1(d)))

	j = sphericaldegree(d)
	v = Jy_eigen!(j, A)
	λ = -j:j
	updatebeta!(d,β)
	djmatrix_fill!(d,j,β,λ,v)
end

function WignerdMatrix(T::Type, j, β::Real)
	d = WignerdMatrix{T}(undef, j, β)
	
	A = Matrix{ComplexF64}(undef, twojp1(j), twojp1(j))
	
	WignerdMatrix!(d,β,A)
end
WignerdMatrix(j, β::Real) = WignerdMatrix(Float64, j, β)

function WignerDMatrix(::Type{Complex{R}}, j, α::Real, β::Real, γ::Real) where {R<:Real}
	d = WignerdMatrix(R, j, β)
	WignerDMatrix(Complex{R}, α, d, γ)
end
WignerDMatrix(j, α::Real, β::Real, γ::Real) = WignerDMatrix(ComplexF64, j, α, β, γ)