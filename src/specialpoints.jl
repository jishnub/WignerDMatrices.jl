abstract type SpecialBetas <: Real end
struct BetaPiby2 <: SpecialBetas end
abstract type ScaledPi <: SpecialBetas end
struct BetaZero <: ScaledPi end
struct BetaPi <: ScaledPi end

# Returns exp(i α π/2) = cos(α π/2) + i*sin(α π/2)
function cis_special(α::HalfInt,::BetaPiby2)
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
		res = cis(α*float(BetaPiby2()))
	end
	return res
end

cis_special(α::HalfInt,::BetaZero) = one(ComplexF64)

# Returns exp(i α π) = cos(α π) + i*sin(α π)
# If α is an integer then this is equal to cos(απ) = (-1)^α
# If α is a half-integer, then this is equal to im*sin(απ) = im*(-1)^(2α)
function cis_special(α::HalfInt,::BetaPi)
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

cis_special(α, p::SpecialBetas) = cis_special(HalfInt(α),p)

Base.cos(::BetaZero) = one(Float64)
Base.cos(::BetaPi) = -one(Float64)
Base.cos(::BetaPiby2) = zero(Float64)
Base.sin(::BetaPiby2) = one(Float64)
Base.sin(::ScaledPi) = zero(Float64)

Base.Float64(::BetaPiby2) = π/2
Base.Float64(::BetaZero) = zero(Float64)
Base.Float64(::BetaPi) = Float64(π)
Base.AbstractFloat(p::SpecialBetas) = Float64(p)

Base.one(::SpecialBetas) = one(Float64)
Base.zero(::SpecialBetas) = zero(Float64)

Base.promote_rule(::Type{<:SpecialBetas},::Type{Float64}) = Float64
Base.promote_rule(::Type{<:SpecialBetas},T::Type{<:Real}) = promote_type(Float64,T)