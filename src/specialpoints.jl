abstract type SpecialAngles <: Real end
struct Piby2Radians <: SpecialAngles end
abstract type ScaledPi <: SpecialAngles end
struct ZeroRadians <: ScaledPi end
struct PiRadians <: ScaledPi end

# Returns exp(i α π/2) = cos(α π/2) + i*sin(α π/2)
function cis_special(α::HalfInt,::Piby2Radians)
	if isinteger(α)
		res = zero(ComplexF64)
		αmod4 = mod(Integer(α),4)
		if αmod4 == 0
			res += one(res)
		elseif αmod4 == 1
			res += one(res)*im
		elseif αmod4 == 2
			res -= one(res)
		elseif αmod4 == 3
			res -= one(res)*im
		end
	else
		res = cis(α*float(Piby2Radians()))
	end
	return res
end

cis_special(α::HalfInt,::ZeroRadians) = one(ComplexF64)

# Returns exp(i α π) = cos(α π) + i*sin(α π)
# If α is an integer then this is equal to cos(απ) = (-1)^α
# If α is a half-integer, then this is equal to im*sin(απ) = im*(-1)^(2α)
function cis_special(α::HalfInt,::PiRadians)
	if isinteger(α)
		res = one(ComplexF64)
		if isodd(Integer(α))
			res *= -1
		end
	else
		res = Complex{Float64}(0,1)
		h = twice(Integer,α)
		if mod(h,4) == 3
			res *= -1
		end
	end
	return res
end

cis_special(α, p::SpecialAngles) = cis_special(HalfInt(α),p)

Base.cos(::ZeroRadians) = one(Float64)
Base.cos(::PiRadians) = -one(Float64)
Base.cos(::Piby2Radians) = zero(Float64)
Base.sin(::Piby2Radians) = one(Float64)
Base.sin(::ScaledPi) = zero(Float64)

Base.Float64(::Piby2Radians) = π/2
Base.Float64(::ZeroRadians) = zero(Float64)
Base.Float64(::PiRadians) = Float64(π)
Base.AbstractFloat(p::SpecialAngles) = Float64(p)

Base.one(::SpecialAngles) = one(Float64)
Base.zero(::SpecialAngles) = zero(Float64)

Base.promote_rule(::Type{<:SpecialAngles},::Type{Float64}) = Float64
Base.promote_rule(::Type{<:SpecialAngles},T::Type{<:Real}) = promote_type(Float64,T)