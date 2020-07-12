abstract type SpecialAngle{T<:Real} <: Real end
struct Piby2Radians{T<:Real} <: SpecialAngle{T} end
abstract type ScaledPi{T<:Real} <: SpecialAngle{T} end
struct ZeroRadians{T<:Real} <: ScaledPi{T} end
struct TwoNPiRadians{T} <: ScaledPi{T}
	n :: Int

	function TwoNPiRadians{T}(n::Int) where {T}
		new{T}(n)
	end
end
struct TwoNPlusOnePiRadians{T} <: ScaledPi{T}
	n :: Int

	function TwoNPlusOnePiRadians{T}(n::Int) where {T}
		new{T}(n)
	end
end

for DT in [:ZeroRadians, :Piby2Radians]
	@eval $DT() = $DT{Float64}()
	@eval $DT{T}(::$DT) where {T} = $DT{T}()
end

for DT in [:TwoNPiRadians, :TwoNPlusOnePiRadians]
	@eval $DT(n::Int) = $DT{Float64}(n)
	@eval $DT{T}(x::$DT) where {T} = $DT{T}(x.n)
end

const Zero = ZeroRadians()
const Pi = TwoNPlusOnePiRadians(0)
const TwoPi = TwoNPiRadians(1)
const ThreePi = TwoNPlusOnePiRadians(1)
const FourPi = TwoNPiRadians(2)
const Piby2 = Piby2Radians()

Base.:(-)(a::SpecialAngle) = -float(a)
Base.:(-)(a::ZeroRadians) = a
Base.:(-)(a::T) where {T<:TwoNPiRadians} = T(-a.n)
Base.:(-)(a::T) where {T<:TwoNPlusOnePiRadians} = T(-a.n - 1)

function Base.:(+)(::ZeroRadians{T1}, ::ZeroRadians{T2}) where {T1,T2}
	ZeroRadians{promote_type(T1,T2)}()
end
for DT in [:Piby2Radians]
	@eval function Base.:(+)(::ZeroRadians{T1}, ::$DT{T2}) where {T1,T2}
		$DT{promote_type(T1,T2)}()
	end
	@eval function Base.:(+)(::$DT{T2}, ::ZeroRadians{T1}) where {T1,T2}
		$DT{promote_type(T1,T2)}()
	end
end

for DT in [:TwoNPiRadians, :TwoNPlusOnePiRadians]
	@eval function Base.:(+)(::ZeroRadians{T1}, x::$DT{T2}) where {T1,T2}
		$DT{promote_type(T1,T2)}(x.n)
	end
end

function Base.:(+)(::Piby2Radians{T1}, ::Piby2Radians{T2}) where {T1,T2}
	TwoNPlusOnePiRadians{promote_type(T1,T2)}(0)
end
function Base.:(+)(a::TwoNPiRadians{T1}, b::TwoNPiRadians{T2}) where {T1,T2}
	TwoNPiRadians{promote_type(T1,T2)}(a.n + b.n)
end

function Base.:(+)(a::TwoNPiRadians{T1}, b::TwoNPlusOnePiRadians{T2}) where {T1,T2}
	TwoNPlusOnePiRadians{promote_type(T1,T2)}(a.n + b.n)
end
function Base.:(+)(a::TwoNPlusOnePiRadians{T2}, b::TwoNPiRadians{T1}) where {T1,T2}
	TwoNPlusOnePiRadians{promote_type(T1,T2)}(a.n + b.n)
end
function Base.:(+)(a::TwoNPlusOnePiRadians{T2}, b::TwoNPlusOnePiRadians{T1}) where {T1,T2}
	TwoNPiRadians{promote_type(T1,T2)}(a.n + b.n + 1)
end

# Returns exp(i α π/2) = cos(α π/2) + i*sin(α π/2)
# if α is a half-integer, we set α = m/2 for an odd integer m
# exp(i m π/4) = cos(m π/4) + i*sin(m π/4)
function cis_special(::Type{T}, α::HalfInt, p::Piby2Radians) where {T}
	if isinteger(α)
		αmod4 = mod(floor(Int,α),4)
		if αmod4 == 0
			res = Complex(one(T),zero(T))
		elseif αmod4 == 1
			res = Complex(zero(T),one(T))
		elseif αmod4 == 2
			res = Complex(-one(T),zero(T))
		elseif αmod4 == 3
			res = Complex(zero(T),-one(T))
		end
	else
		invsqrt2 = 1/√(T(2))
		m = mod(floor(Int,twice(α)),8)
		if m == 1
			res = Complex(invsqrt2, invsqrt2)
		elseif m == 3
			res = Complex(-invsqrt2, invsqrt2)
		elseif m == 5
			res = Complex(-invsqrt2, -invsqrt2)
		elseif m == 7
			res = Complex(invsqrt2, -invsqrt2)
		end
	end
	return res
end

cis_special(::Type{T}, α::HalfInt, ::ZeroRadians) where {T} = one(Complex{T})

# Returns exp(i α 2πn) = cos(2πn α) + i*sin(2πn α) = (-1)^(n * 2α)
function cis_special(::Type{T}, α::HalfInt, x::TwoNPiRadians) where {T}
	isodd(x.n * twice(α)) ? Complex(-one(T)) : Complex(one(T))
end

# Returns exp(i α (2n+1)π) = cos((2n+1)π α) + i*sin((2n+1)π α)
# if α is an integer, we obtain exp(i α (2n+1)π) = cos((2n+1)π α) = (-1)^((2n+1)α)
# if α is a half-integer = m/2 for odd m, we obtain 
# exp(i α (2n+1)π) = i*sin((2n+1) * m * π/2)
function cis_special(::Type{T}, α::HalfInt, x::TwoNPlusOnePiRadians) where {T}
	r, i = zero(T), zero(T)
	if isinteger(α)
		if isodd(floor(Int,α))
			r = -one(T)
		else
			r = one(T)
		end
	else
		m = twice(α)
		p = (2x.n + 1)*m
		if mod(p, 4) == 1
			i = one(T)
		else
			i = -one(T)
		end
	end
	Complex{T}(r, i)
end

cis_special(α::Real, p::SpecialAngle{T}) where {T} = cis_special(T, HalfInt(α), p)

Base.cos(::ZeroRadians{T}) where {T} = one(T)
Base.cos(::TwoNPiRadians{T}) where {T} = one(T)
Base.cos(::TwoNPlusOnePiRadians{T}) where {T} = -one(T)
Base.cos(::Piby2Radians{T}) where {T} = zero(T)
Base.sin(::Piby2Radians{T}) where {T} = one(T)
Base.sin(::ScaledPi{T}) where {T} = zero(T)

for T in [:Float16, :Float32, :Float64, :BigFloat]
	@eval Base.$T(::Piby2Radians) = $T(π)/2
	@eval Base.$T(::ZeroRadians) = zero($T)
	@eval Base.$T(x::TwoNPiRadians) = $T(π)*2x.n
	@eval Base.$T(x::TwoNPlusOnePiRadians) = $T(π)*(2x.n + 1)
end
Base.AbstractFloat(p::SpecialAngle{T}) where {T} = T(p)

Base.one(::SpecialAngle{T}) where {T} = one(T)
Base.zero(::SpecialAngle{T}) where {T} = zero(T)

Base.promote_rule(::Type{<:SpecialAngle{T}}, R::Type{<:Real}) where {T} = promote_type(T,R)

Base.show(io::IO, x::ZeroRadians) = print(io, "0")
Base.show(io::IO, x::Piby2Radians) = print(io, "π/2")
function Base.show(io::IO, x::TwoNPiRadians)
	if iszero(x.n)
		print(io, 0)
	else
		print(io, 2x.n,"π")
	end
end
function Base.show(io::IO, x::TwoNPlusOnePiRadians)
	if iszero(x.n)
		print(io, "π")
	elseif x.n == -1
		print(io, "-π")
	else
		print(io, 2x.n + 1,"π")
	end
end