using WignerDMatrices
using HalfIntegers
using LinearAlgebra
using Test

import WignerDMatrices: cis_special, Piby2Radians, ZeroRadians,
TwoNPiRadians, TwoNPlusOnePiRadians
import WignerDMatrices: WignerdMatrixContainer, flatind, filledelements, 
flatind_phase, sphericaldegree

@test isempty(Test.detect_ambiguities(Base, Core, WignerDMatrices))

@testset "special points" begin
	@testset "float" begin
		@testset "Piby2Radians" begin
		    @test float(Piby2Radians()) == π/2
	    	@test AbstractFloat(Piby2Radians()) == π/2
	    	@test Float64(Piby2Radians()) == π/2
	    	@test float(Piby2Radians{BigFloat}()) == big(π)/2
	    	@test Float64(Piby2Radians{BigFloat}()) == π/2
		end
		@testset "ZeroRadians" begin
			@test float(ZeroRadians()) == 0
	    	@test AbstractFloat(ZeroRadians()) == 0
	    	@test Float64(ZeroRadians()) == 0
	    	@test float(ZeroRadians{BigFloat}()) isa BigFloat
	    	@test float(ZeroRadians{BigFloat}()) == big(0.0)
		end
		@testset "TwoNPiRadians" begin
		    @test float(TwoNPiRadians(1)) == 2π
		    @test float(TwoNPiRadians{BigFloat}(1)) == big(2)*π
		    @test float(TwoNPiRadians(-1)) == -2π
		    @test float(TwoNPiRadians{BigFloat}(-1)) == big(-2)*π
		end
		@testset "TwoNPlusOnePiRadians" begin
		   	@test float(TwoNPlusOnePiRadians(1)) == 3π
		    @test float(TwoNPlusOnePiRadians{BigFloat}(1)) == big(3)*π
		   	@test float(TwoNPlusOnePiRadians(-1)) == -π
		    @test float(TwoNPlusOnePiRadians{BigFloat}(-1)) == big(-1)*π
		end
		@testset "TwoNPiRadians and ZeroRadians" begin
		    @test float(TwoNPiRadians(0)) == float(ZeroRadians())
		    @test float(TwoNPiRadians{BigFloat}(0)) == float(ZeroRadians{BigFloat}())
		end
	end
	@testset "promote" begin
		for T in [Int,Float32,Float64]
		    @test promote_rule(Piby2Radians{Float64},T) == Float64
		    @test promote_rule(ZeroRadians{Float64},T) == Float64
		end
		for T in [Int,Float32,Float64]
		    @test promote_rule(Piby2Radians{BigFloat},T) == BigFloat
		    @test promote_rule(ZeroRadians{BigFloat},T) == BigFloat
		end
		for T in [BigFloat]
			@test promote_rule(Piby2Radians{Float64},T) == BigFloat
			@test promote_rule(ZeroRadians{Float64},T) == BigFloat
		end
	end
	@testset "zero and one" begin
		function test(T, n...)
			@test zero(T(n...)) == zero(Float64)
	    	@test one(T(n...)) == one(Float64)
		    @test zero(T{BigFloat}(n...)) isa BigFloat
		    @test zero(T{BigFloat}(n...)) == zero(BigFloat)
	    	@test one(T{BigFloat}(n...)) isa BigFloat
	    	@test one(T{BigFloat}(n...)) == one(BigFloat)
		end
		@testset "Piby2Radians" begin
		    test(Piby2Radians)
		end
		@testset "ZeroRadians" begin
			test(ZeroRadians)
		end
		@testset "TwoNPiRadians" begin
		   	test(TwoNPiRadians, 1)
		end
		@testset "TwoNPlusOnePiRadians" begin
		   	test(TwoNPlusOnePiRadians, 1) 
		end
	end
	@testset "negate" begin
	    @test -ZeroRadians() == ZeroRadians()
	    @test -Piby2Radians() == -π/2
	    @test -Piby2Radians{BigFloat}() == -big(π)/2
	    @test -TwoNPlusOnePiRadians{BigFloat}(0) == -big(π)
	    @test -TwoNPiRadians(1) == -2π
	    @test -TwoNPiRadians{BigFloat}(1) == -2big(π)
	    @test -TwoNPlusOnePiRadians(1) == -3π
	    @test -TwoNPlusOnePiRadians{BigFloat}(1) == -3big(π)
	end
	@testset "add" begin
		@testset "ZeroRadians" begin
			@test ZeroRadians() + ZeroRadians() == ZeroRadians()
	    	@test ZeroRadians{Float64}() + ZeroRadians{BigFloat}() === ZeroRadians{BigFloat}()
	    	@test ZeroRadians{BigFloat}() + ZeroRadians{Float64}() === ZeroRadians{BigFloat}()

	    	@test ZeroRadians() + Piby2Radians() == Piby2Radians()
	    	@test ZeroRadians() + Piby2Radians{BigFloat}() == Piby2Radians{BigFloat}()
	    	@test ZeroRadians{BigFloat}() + Piby2Radians() == Piby2Radians{BigFloat}()
	    	
	    	@test ZeroRadians() + TwoNPiRadians(1) == TwoNPiRadians(1)
	    	@test ZeroRadians() + TwoNPiRadians{BigFloat}(1) == TwoNPiRadians{BigFloat}(1)
	    	@test ZeroRadians{BigFloat}() + TwoNPiRadians(1) == TwoNPiRadians{BigFloat}(1)
	    	
	    	@test ZeroRadians() + TwoNPlusOnePiRadians(1) == TwoNPlusOnePiRadians(1)
	    	@test ZeroRadians() + TwoNPlusOnePiRadians{BigFloat}(1) == TwoNPlusOnePiRadians{BigFloat}(1)
	    	@test ZeroRadians{BigFloat}() + TwoNPlusOnePiRadians(1) == TwoNPlusOnePiRadians{BigFloat}(1)
		end
		@testset "Piby2Radians" begin
			@test Piby2Radians() + Piby2Radians() == TwoNPlusOnePiRadians(0)
	    	@test Piby2Radians{Float64}() + Piby2Radians{BigFloat}() === TwoNPlusOnePiRadians{BigFloat}(0)
	    	@test Piby2Radians{BigFloat}() + Piby2Radians{Float64}() === TwoNPlusOnePiRadians{BigFloat}(0)

	    	@test Piby2Radians() + ZeroRadians() == Piby2Radians()
	    	@test Piby2Radians{BigFloat}() + ZeroRadians() == Piby2Radians{BigFloat}()
	    	@test Piby2Radians() + ZeroRadians{BigFloat}() == Piby2Radians{BigFloat}()

	    	@test Piby2Radians() + TwoNPiRadians(1) == 5π/2
	    	@test Piby2Radians() + TwoNPiRadians{BigFloat}(1) == 5big(π)/2
	    	@test Piby2Radians{BigFloat}() + TwoNPiRadians(1) == 5big(π)/2

	    	@test Piby2Radians() + TwoNPlusOnePiRadians(1) == 7π/2
	    	@test Piby2Radians{BigFloat}() + TwoNPlusOnePiRadians(1) == 7big(π)/2
	    	@test Piby2Radians() + TwoNPlusOnePiRadians{BigFloat}(1) == 7big(π)/2
		end
		@testset "TwoNPiRadians" begin
		    @test TwoNPiRadians(1) + TwoNPiRadians(1) == TwoNPiRadians(2)
	    	@test TwoNPiRadians{Float64}(1) + TwoNPiRadians{BigFloat}(1) === TwoNPiRadians{BigFloat}(2)
	    	@test TwoNPiRadians{BigFloat}(1) + TwoNPiRadians{Float64}(1) === TwoNPiRadians{BigFloat}(2)

	    	@test TwoNPiRadians(1) + ZeroRadians() == TwoNPiRadians(1)
	    	@test TwoNPiRadians{BigFloat}(1) + ZeroRadians() == TwoNPiRadians{BigFloat}(1)
	    	@test TwoNPiRadians(1) + ZeroRadians{BigFloat}() == TwoNPiRadians{BigFloat}(1)

	    	@test TwoNPiRadians(1) + Piby2Radians() == 5π/2
	    	@test TwoNPiRadians{BigFloat}(1) + Piby2Radians() == 5big(π)/2
	    	@test TwoNPiRadians(1) + Piby2Radians{BigFloat}() == 5big(π)/2

	    	@test TwoNPiRadians(1) + TwoNPlusOnePiRadians(0) == TwoNPlusOnePiRadians(1)
	    	@test TwoNPiRadians{BigFloat}(1) + TwoNPlusOnePiRadians(0) == TwoNPlusOnePiRadians{BigFloat}(1)
	    	@test TwoNPiRadians(1) + TwoNPlusOnePiRadians{BigFloat}(0) == TwoNPlusOnePiRadians{BigFloat}(1)

	    	@test TwoNPiRadians(1) + TwoNPlusOnePiRadians(1) == TwoNPlusOnePiRadians(2)
	    	@test TwoNPiRadians{BigFloat}(1) + TwoNPlusOnePiRadians(1) == TwoNPlusOnePiRadians{BigFloat}(2)
	    	@test TwoNPiRadians(1) + TwoNPlusOnePiRadians{BigFloat}(1) == TwoNPlusOnePiRadians{BigFloat}(2)
		end
		@testset "TwoNPlusOnePiRadians" begin
		    @test TwoNPlusOnePiRadians(1) + TwoNPlusOnePiRadians(1) == TwoNPiRadians(3)
	    	@test TwoNPlusOnePiRadians{Float64}(1) + TwoNPlusOnePiRadians{BigFloat}(1) === TwoNPiRadians{BigFloat}(3)
	    	@test TwoNPlusOnePiRadians{BigFloat}(1) + TwoNPlusOnePiRadians{Float64}(1) === TwoNPiRadians{BigFloat}(3)

	    	@test TwoNPlusOnePiRadians(1) + ZeroRadians() == TwoNPlusOnePiRadians(1)
	    	@test TwoNPlusOnePiRadians{BigFloat}(1) + ZeroRadians() == TwoNPlusOnePiRadians{BigFloat}(1)
	    	@test TwoNPlusOnePiRadians(1) + ZeroRadians{BigFloat}() == TwoNPlusOnePiRadians{BigFloat}(1)

	    	@test TwoNPlusOnePiRadians(1) + Piby2Radians() == 7π/2
	    	@test TwoNPlusOnePiRadians{BigFloat}(1) + Piby2Radians() == 7big(π)/2
	    	@test TwoNPlusOnePiRadians(1) + Piby2Radians{BigFloat}() == 7big(π)/2

	    	@test TwoNPlusOnePiRadians(1) + TwoNPlusOnePiRadians(0) == TwoNPiRadians(2)
	    	@test TwoNPlusOnePiRadians{BigFloat}(1) + TwoNPlusOnePiRadians(0) == TwoNPiRadians{BigFloat}(2)
	    	@test TwoNPlusOnePiRadians(1) + TwoNPlusOnePiRadians{BigFloat}(0) == TwoNPiRadians{BigFloat}(2)

	    	@test TwoNPlusOnePiRadians(1) + TwoNPiRadians(1) == TwoNPlusOnePiRadians(2)
	    	@test TwoNPlusOnePiRadians(1) + TwoNPiRadians{BigFloat}(1) == TwoNPlusOnePiRadians{BigFloat}(2)
	    	@test TwoNPlusOnePiRadians{BigFloat}(1) + TwoNPiRadians(1) == TwoNPlusOnePiRadians{BigFloat}(2)
		end
	end
	@testset "trigonometric functions" begin
	    @testset "Piby2Radians" begin
		    for α = -10:1//2:10
		    	@test cis_special(HalfInt(α),Piby2Radians()) isa ComplexF64
		    	@test cis_special(HalfInt(α),Piby2Radians()) ≈ cis(α*π/2)
		    	@test cis_special(HalfInt(α),Piby2Radians{BigFloat}()) isa Complex{BigFloat}
		    	@test cis_special(HalfInt(α),Piby2Radians{BigFloat}()) ≈ cis(α*big(π)/2)
		    end
	    	@test cos(Piby2Radians()) == 0
	    	@test cos(Piby2Radians{BigFloat}()) isa BigFloat
	    	@test cos(Piby2Radians{BigFloat}()) == 0
	    	@test sin(Piby2Radians()) == 1
	    	@test sin(Piby2Radians{BigFloat}()) isa BigFloat
	    	@test sin(Piby2Radians{BigFloat}()) == 1
		end
		@testset "ZeroRadians" begin
			for α = -10:1//2:10
		    	@test cis_special(HalfInt(α),ZeroRadians()) isa ComplexF64
		    	@test cis_special(HalfInt(α),ZeroRadians()) ≈ cis(0)
		    	@test cis_special(HalfInt(α),ZeroRadians{BigFloat}()) isa Complex{BigFloat}
		    	@test cis_special(HalfInt(α),ZeroRadians{BigFloat}()) ≈ cis(big(0))
		    end
		    @test cos(ZeroRadians()) == 1
		    @test cos(ZeroRadians{BigFloat}()) isa BigFloat
		    @test cos(ZeroRadians{BigFloat}()) == 1
	    	@test sin(ZeroRadians()) == 0
	    	@test sin(ZeroRadians{BigFloat}()) isa BigFloat
	    	@test sin(ZeroRadians{BigFloat}()) == 0
		end
		@testset "TwoNPiRadians" begin
		    for α = -10:1//2:10
		    	@test cis_special(HalfInt(α),TwoNPiRadians(1)) isa ComplexF64
		    	@test cis_special(HalfInt(α),TwoNPiRadians(1)) ≈ cis(α*2π)
		    	@test cis_special(HalfInt(α),TwoNPiRadians{BigFloat}(1)) isa Complex{BigFloat}
		    	@test cis_special(HalfInt(α),TwoNPiRadians{BigFloat}(1)) ≈ cis(α*2big(π))
		    end
		    for n = -3:3
			    @test cos(TwoNPiRadians(n)) == 1
			    @test cos(TwoNPiRadians{BigFloat}(n)) isa BigFloat
			    @test cos(TwoNPiRadians{BigFloat}(n)) == 1
		    	@test sin(TwoNPiRadians(n)) == 0
		    	@test sin(TwoNPiRadians{BigFloat}(n)) isa BigFloat
		    	@test sin(TwoNPiRadians{BigFloat}(n)) == 0
		    end
		end
		@testset "TwoNPlusOnePiRadians" begin
		    for α = -10:1//2:10
		    	@test cis_special(HalfInt(α),TwoNPlusOnePiRadians(1)) isa ComplexF64
		    	@test cis_special(HalfInt(α),TwoNPlusOnePiRadians(1)) ≈ cis(α*3π)
		    	@test cis_special(HalfInt(α),TwoNPlusOnePiRadians{BigFloat}(1)) isa Complex{BigFloat}
		    	@test cis_special(HalfInt(α),TwoNPlusOnePiRadians{BigFloat}(1)) ≈ cis(α*3big(π))
		    end
		    for n = -3:3
			    @test cos(TwoNPlusOnePiRadians(n)) == -1
			    @test cos(TwoNPlusOnePiRadians{BigFloat}(n)) isa BigFloat
			    @test cos(TwoNPlusOnePiRadians{BigFloat}(n)) == -1
		    	@test sin(TwoNPlusOnePiRadians(n)) == 0
		    	@test sin(TwoNPlusOnePiRadians{BigFloat}(n)) isa BigFloat
		    	@test sin(TwoNPlusOnePiRadians{BigFloat}(n)) == 0
		    end
		end
	end
	@testset "wignerdmatrixelement" begin
		function testapprox(m,n,dj_m_n,dj_m_n2)
	    	@test begin 
	    		res = isapprox(dj_m_n,dj_m_n2,atol=1e-13,rtol=sqrt(eps(Float64)))
	    		if !res
	    			@show m n dj_m_n dj_m_n2
	    		end
	    		res
	    	end
	    end
	    function tests(j,v)
	    	@testset "Piby2Radians" begin
		        for m in -j:j, n in -j:j
		        	dj_m_n = WignerDMatrices.wignerdmatrixelement(j,m,n,π/2 ,v)
		        	dj_m_n2 = WignerDMatrices.wignerdmatrixelement(j,m,n,Piby2Radians(),v)

		        	testapprox(m,n,dj_m_n,dj_m_n2)
		        end
		    end
		    @testset "ZeroRadians" begin
		        for m in -j:j, n in -j:j
		        	dj_m_n = WignerDMatrices.wignerdmatrixelement(j,m,n,0,v)
		        	dj_m_n2 = WignerDMatrices.wignerdmatrixelement(j,m,n,ZeroRadians(),v)

		        	testapprox(m,n,dj_m_n,dj_m_n2)
		        end
		    end
		    @testset "TwoNPiRadians" begin
		        for m in -j:j, n in -j:j
		        	dj_m_n = WignerDMatrices.wignerdmatrixelement(j,m,n,2π,v)
		        	dj_m_n2 = WignerDMatrices.wignerdmatrixelement(j,m,n,TwoNPiRadians(1),v)

		        	testapprox(m,n,dj_m_n,dj_m_n2)

		        	dj_m_n = WignerDMatrices.wignerdmatrixelement(j,m,n,0,v)
		        	dj_m_n2 = WignerDMatrices.wignerdmatrixelement(j,m,n,TwoNPiRadians(0),v)

		        	testapprox(m,n,dj_m_n,dj_m_n2)
		        	
		        	dj_m_n = WignerDMatrices.wignerdmatrixelement(j,m,n,ZeroRadians(),v)
		        	dj_m_n2 = WignerDMatrices.wignerdmatrixelement(j,m,n,TwoNPiRadians(0),v)

		        	testapprox(m,n,dj_m_n,dj_m_n2)
		        end
		    end
		    @testset "TwoNPlusOnePiRadians" begin
		        for m in -j:j, n in -j:j
		        	dj_m_n = WignerDMatrices.wignerdmatrixelement(j,m,n,3π,v)
		        	dj_m_n2 = WignerDMatrices.wignerdmatrixelement(j,m,n,WignerDMatrices.ThreePi,v)

		        	dj_m_n = WignerDMatrices.wignerdmatrixelement(j,m,n,π,v)
		        	dj_m_n2 = WignerDMatrices.wignerdmatrixelement(j,m,n,WignerDMatrices.Pi,v)

		        	testapprox(m,n,dj_m_n,dj_m_n2)
		        end
		    end
		end

		@testset "integer j" begin
		    for j = 0:10
				A = zeros(ComplexF64,2j+1,2j+1)
			    v = WignerDMatrices.eigenvecsJy!(j,A)

			    tests(j,v)
			end
		end
		@testset "half-integer j" begin
			for j = 1//2:19//2
				A = zeros(ComplexF64,Int(2j+1),Int(2j+1))
			    v = WignerDMatrices.eigenvecsJy!(HalfInt(j),A)

		    	tests(j,v)
		    end
		end
	end
end
@testset "WignerdMatrixContainer" begin
    @testset "flatind" begin
    	@testset "j = 1/2" begin
    	    j = HalfInt(1/2)
    	    @test filledelements(j) == 2
    	    @test flatind(j, 1, 1) == 1
            @test flatind(j, 2, 1) == 2
    	end
    	@testset "j = 1" begin
        	j = 1
        	@test filledelements(j) == 4
            @test flatind(j, 1, 1) == 1
            @test flatind(j, 2, 1) == 2
            @test flatind(j, 3, 1) == 3
            @test flatind(j, 2, 2) == 4
        end
    	@testset "j = 3/2" begin
        	j = HalfInt(3/2)
        	@test filledelements(j) == 6
            @test flatind(j, 1, 1) == 1
            @test flatind(j, 2, 1) == 2
            @test flatind(j, 3, 1) == 3
            @test flatind(j, 4, 1) == 4
            @test flatind(j, 2, 2) == 5
            @test flatind(j, 3, 2) == 6
        end
        @testset "j = 2" begin
        	j = 2
        	@test filledelements(j) == 9
            @test flatind(j, 1, 1) == 1
            @test flatind(j, 2, 1) == 2
            @test flatind(j, 3, 1) == 3
            @test flatind(j, 4, 1) == 4
            @test flatind(j, 5, 1) == 5
            @test flatind(j, 2, 2) == 6
            @test flatind(j, 3, 2) == 7
            @test flatind(j, 4, 2) == 8
            @test flatind(j, 3, 3) == 9
        end
    end
    @testset "flating phase" begin
        @testset "j = 1/2" begin
			j = HalfInt(1/2)
			@testset "saved" begin
				@test flatind_phase(j, 1, 1) == (1, 1)
				@test flatind_phase(j, 2, 1) == (2, 1)
			end
			@testset "right" begin
				@test flatind_phase(j, 1, 2) == (2, -1)
				@test flatind_phase(j, 2, 2) == (1, 1)
			end
    	end
    	@testset "j = 1" begin
        	j = 1
        	@testset "saved" begin
	            @test flatind_phase(j, 1, 1) == (1, 1)
	            @test flatind_phase(j, 2, 1) == (2, 1)
	            @test flatind_phase(j, 3, 1) == (3, 1)
            	@test flatind_phase(j, 2, 2) == (4, 1)
        	end
        	@testset "top" begin
            	@test flatind_phase(j, 1, 2) == (2, -1)
        	end
        	@testset "bottom" begin
            	@test flatind_phase(j, 3, 2) == (2, 1)
        	end
        	@testset "right" begin
            	@test flatind_phase(j, 1, 3) == (3, 1)
	            @test flatind_phase(j, 2, 3) == (2, -1)
	            @test flatind_phase(j, 3, 3) == (1, 1)
        	end
        end
    	@testset "j = 3/2" begin
        	j = HalfInt(3/2)
        	@testset "saved" begin
	        	@test flatind_phase(j, 1, 1) == (1, 1)
	            @test flatind_phase(j, 2, 1) == (2, 1)
	            @test flatind_phase(j, 3, 1) == (3, 1)
	            @test flatind_phase(j, 4, 1) == (4, 1)
	            @test flatind_phase(j, 2, 2) == (5, 1)
	            @test flatind_phase(j, 3, 2) == (6, 1)
        	end
        	@testset "top" begin
            	@test flatind_phase(j, 1, 2) == (2, -1)
            	@test flatind_phase(j, 1, 3) == (3, 1)
        	end
        	@testset "bottom" begin
            	@test flatind_phase(j, 4, 2) == (3, 1)
            	@test flatind_phase(j, 4, 3) == (2, 1)
        	end
        	@testset "right" begin
	            @test flatind_phase(j, 2, 3) == (6, -1)
	            @test flatind_phase(j, 3, 3) == (5, 1)
	            @test flatind_phase(j, 1, 4) == (4, -1)
	            @test flatind_phase(j, 2, 4) == (3, 1)
	            @test flatind_phase(j, 3, 4) == (2, -1)
	            @test flatind_phase(j, 4, 4) == (1, 1)
        	end
        end
        @testset "j = 2" begin
        	j = 2
        	@testset "saved" begin
	        	@test flatind_phase(j, 1, 1) == (1, 1)
	            @test flatind_phase(j, 2, 1) == (2, 1)
	            @test flatind_phase(j, 3, 1) == (3, 1)
	            @test flatind_phase(j, 4, 1) == (4, 1)
	            @test flatind_phase(j, 5, 1) == (5, 1)
	            @test flatind_phase(j, 2, 2) == (6, 1)
	            @test flatind_phase(j, 3, 2) == (7, 1)
	            @test flatind_phase(j, 4, 2) == (8, 1)
            	@test flatind_phase(j, 3, 3) == (9, 1)
        	end
        	@testset "top" begin
            	@test flatind_phase(j, 1, 2) == (2, -1)
            	@test flatind_phase(j, 1, 3) == (3, 1)
            	@test flatind_phase(j, 1, 4) == (4, -1)
            	@test flatind_phase(j, 2, 3) == (7, -1)
        	end
        	@testset "bottom" begin
            	@test flatind_phase(j, 5, 2) == (4, 1)
            	@test flatind_phase(j, 5, 3) == (3, 1)
            	@test flatind_phase(j, 5, 4) == (2, 1)
            	@test flatind_phase(j, 4, 3) == (7, 1)
        	end
        	@testset "right" begin
	            @test flatind_phase(j, 2, 4) == (8, 1)
	            @test flatind_phase(j, 3, 4) == (7, -1)
	            @test flatind_phase(j, 4, 4) == (6, 1)
	            @test flatind_phase(j, 1, 5) == (5, 1)
	            @test flatind_phase(j, 2, 5) == (4, -1)
	            @test flatind_phase(j, 3, 5) == (3, 1)
	            @test flatind_phase(j, 4, 5) == (2, -1)
	            @test flatind_phase(j, 5, 5) == (1, 1)
        	end
        end
    end
end
@testset "Jy" begin
	jmax = 2
    A = zeros(ComplexF64,2jmax+1,2jmax+1)

    @testset "Jy_matrix_zbasis!" begin
    	j = 0
	    h = WignerDMatrices.Jy_matrix_zbasis!(j,A)
	    @test size(h) == (2j+1,2j+1)
	    @test ishermitian(h)
	    @test iszero(h[1])

	    j = 1
	    h = WignerDMatrices.Jy_matrix_zbasis!(j,A)
	    @test size(h) == (2j+1,2j+1)
	    @test ishermitian(h)
	    @test all(iszero,diag(h))
	    @test h[1,2] == im*WignerDMatrices.X(j,-j+1)/2
	    @test h[2,3] == im*WignerDMatrices.X(j,-j+2)/2
        
        j = 2
        h = WignerDMatrices.Jy_matrix_zbasis!(j,A)
	    @test size(h) == (2j+1,2j+1)
	    @test ishermitian(h)
	    @test all(iszero,diag(h))
        
        # check that this runs
        v = WignerDMatrices.Jy_eigen_nomatch!(1, A)
    end
end
@testset "WignerdMatrix" begin
	@testset "Constructor" begin
		β = π/3
		@testset "integer j" begin
			j = 2
			jh = HalfInt(j)
			d = WignerdMatrix(j,β)
			@test WignerDMatrices.sphericaldegree(d) == j
			@test d.beta == β
			@test WignerDMatrices.eulerangles(d) == (zero(β),β,zero(β))
			@test axes(d) == (-jh:jh,-jh:jh)
			@test axes(d,1) == axes(d,2) == -jh:jh
			@test size(d) == (Integer(2j+1),Integer(2j+1))
			@test size(d,1) == size(d,2) == Integer(2j+1)

			@test d == d

			for i in eachindex(d.dj)
				@test d[i] == d.dj[i]
			end
		end
		@testset "half-integer j" begin
			j = 3//2
			jh = HalfInt(j)
			d = WignerdMatrix(j,β)
			@test WignerDMatrices.sphericaldegree(d) == j
			@test d.beta == β
			@test WignerDMatrices.eulerangles(d) == (zero(β),β,zero(β))
			@test axes(d) == (-jh:jh,-jh:jh)
			@test axes(d,1) == axes(d,2) == -jh:jh
			@test size(d) == (Integer(2j+1),Integer(2j+1))
			@test size(d,1) == size(d,2) == Integer(2j+1)

			for i in eachindex(d.dj)
				@test d[i] == d.dj[i]
			end
		end
	end
	function checkdjvalues(test,j,n)
		for β in LinRange(0, 4π, n)
			d = WignerdMatrix(j,β)
			test(d,β)
		end

		@testset "Piby2Radians" begin
			β = Piby2Radians()
			d = WignerdMatrix(j,β)
			test(d,β) 
		end
		@testset "ZeroRadians" begin
			β = ZeroRadians()
			d = WignerdMatrix(j,β)
			test(d,β) 
		end
		@testset "TwoNPiRadians" begin
		    for n in -2:2
		   		β = TwoNPiRadians(n)
				d = WignerdMatrix(j,β)
				test(d,β)
			end 	
		end
		@testset "TwoNPlusOnePiRadians" begin
		    for n in -2:2
		   		β = TwoNPlusOnePiRadians(n)
				d = WignerdMatrix(j,β)
				test(d,β)
			end 	
		end
	end
	@testset "d1//2_mn(β)" begin
		function test(d,β)
			@test isapprox(d[-1/2,-1/2],cos(β/2),atol=1e-14,rtol=1e-8)
			@test isapprox(d[1/2,-1/2],-sin(β/2),atol=1e-14,rtol=1e-8)

			@test isapprox(d[-1/2,1/2],sin(β/2),atol=1e-14,rtol=1e-8)
			@test isapprox(d[1/2,1/2],cos(β/2),atol=1e-14,rtol=1e-8)
		end
		
		n = 100
		j = 1//2
		checkdjvalues(test,j,n)
	end
    @testset "d1_mn(β)" begin
		function test(d,β)
			@test isapprox(d[1,1],(1+cos(β))/2,atol=1e-14,rtol=1e-8)
			@test isapprox(d[1,0],-sin(β)/√2,atol=1e-14,rtol=1e-8)
			@test isapprox(d[1,-1],(1-cos(β))/2,atol=1e-14,rtol=1e-8)

			@test isapprox(d[0,1],sin(β)/√2,atol=1e-14,rtol=1e-8)
			@test isapprox(d[0,0],cos(β),atol=1e-14,rtol=1e-8)
			@test isapprox(d[0,-1],-sin(β)/√2,atol=1e-14,rtol=1e-8)

			@test isapprox(d[-1,1],(1-cos(β))/2,atol=1e-14,rtol=1e-8)
			@test isapprox(d[-1,0],sin(β)/√2,atol=1e-14,rtol=1e-8)
			@test isapprox(d[-1,-1],(1+cos(β))/2,atol=1e-14,rtol=1e-8)
		end
		
		n = 100
		j = 1
		checkdjvalues(test,j,n)
	end
	@testset "d3//2_mn(β)" begin
	    function test(d,β)
			@test isapprox(d[-3/2,-3/2],cos(β/2)^3,atol=1e-14,rtol=1e-8)
			@test isapprox(d[-1/2,-3/2],-√3*sin(β/2)*cos(β/2)^2,atol=1e-14,rtol=1e-8)
			@test isapprox(d[1/2,-3/2],√3*sin(β/2)^2*cos(β/2),atol=1e-14,rtol=1e-8)
			@test isapprox(d[3/2,-3/2],-sin(β/2)^3,atol=1e-14,rtol=1e-8)

			@test isapprox(d[-3/2,-1/2],√3*sin(β/2)*cos(β/2)^2,atol=1e-14,rtol=1e-8)
			@test isapprox(d[-1/2,-1/2],cos(β/2)*(3cos(β/2)^2-2),atol=1e-14,rtol=1e-8)
			@test isapprox(d[1/2,-1/2],sin(β/2)*(3sin(β/2)^2-2),atol=1e-14,rtol=1e-8)
			@test isapprox(d[3/2,-1/2],√3*sin(β/2)^2*cos(β/2),atol=1e-14,rtol=1e-8)

			@test isapprox(d[-3/2,1/2],√3*sin(β/2)^2*cos(β/2),atol=1e-14,rtol=1e-8)
			@test isapprox(d[-1/2,1/2],-sin(β/2)*(3sin(β/2)^2-2),atol=1e-14,rtol=1e-8)
			@test isapprox(d[1/2,1/2],cos(β/2)*(3cos(β/2)^2-2),atol=1e-14,rtol=1e-8)
			@test isapprox(d[3/2,1/2],-√3*sin(β/2)*cos(β/2)^2,atol=1e-14,rtol=1e-8)

			@test isapprox(d[-3/2,3/2],sin(β/2)^3,atol=1e-14,rtol=1e-8)
			@test isapprox(d[-1/2,3/2],√3*sin(β/2)^2*cos(β/2),atol=1e-14,rtol=1e-8)
			@test isapprox(d[1/2,3/2],√3*sin(β/2)*cos(β/2)^2,atol=1e-14,rtol=1e-8)
			@test isapprox(d[3/2,3/2],cos(β/2)^3,atol=1e-14,rtol=1e-8)
		end
		
		n = 100
		j = 3//2
		checkdjvalues(test,j,n)
	end
	@testset "d2_mn(β)" begin
		function test(d::WignerdMatrix,β)
			@test isapprox(d[2,1],-sin(β)*(1+cos(β))/2,atol=1e-14,rtol=1e-8)
			@test isapprox(d[2,0],1/2*√(3/2)*sin(β)^2,atol=1e-14,rtol=1e-8)
			@test isapprox(d[2,-1],-sin(β)*(1-cos(β))/2,atol=1e-14,rtol=1e-8)
			
			@test isapprox(d[1,1],(2cos(β)^2+cos(β)-1)/2,atol=1e-14,rtol=1e-8)
			@test isapprox(d[1,0],-√(3/2)*sin(β)*cos(β),atol=1e-14,rtol=1e-8)
			@test isapprox(d[1,-1],-(2cos(β)^2-cos(β)-1)/2,atol=1e-14,rtol=1e-8)

			@test isapprox(d[0,1],√(3/2)*sin(β)*cos(β),atol=1e-14,rtol=1e-8)
			@test isapprox(d[0,0],1/2*(3cos(β)^2-1),atol=1e-14,rtol=1e-8)
			@test isapprox(d[0,-1],-√(3/2)*sin(β)*cos(β),atol=1e-14,rtol=1e-8)

			@test isapprox(d[-1,1],-(2cos(β)^2-cos(β)-1)/2,atol=1e-14,rtol=1e-8)
			@test isapprox(d[-1,0],√(3/2)*sin(β)*cos(β),atol=1e-14,rtol=1e-8)
			@test isapprox(d[-1,-1],(2cos(β)^2+cos(β)-1)/2,atol=1e-14,rtol=1e-8)

			@test isapprox(d[-2,1],sin(β)*(1-cos(β))/2,atol=1e-14,rtol=1e-8)
			@test isapprox(d[-2,0],1/2*√(3/2)*sin(β)^2,atol=1e-14,rtol=1e-8)
			@test isapprox(d[-2,-1],sin(β)*(1+cos(β))/2,atol=1e-14,rtol=1e-8)

			@test isapprox(d[2,2],(1+cos(β))^2/4,atol=1e-14,rtol=1e-8)
			@test isapprox(d[2,-2],(1-cos(β))^2/4,atol=1e-14,rtol=1e-8)
			
			@test isapprox(d[1,2],sin(β)*(1+cos(β))/2,atol=1e-14,rtol=1e-8)
			@test isapprox(d[1,-2],-sin(β)*(1-cos(β))/2,atol=1e-14,rtol=1e-8)
			
			@test isapprox(d[0,2],1/2*√(3/2)*sin(β)^2,atol=1e-14,rtol=1e-8)
			@test isapprox(d[0,-2],1/2*√(3/2)*sin(β)^2,atol=1e-14,rtol=1e-8)

			@test isapprox(d[-1,2],sin(β)*(1-cos(β))/2,atol=1e-14,rtol=1e-8)
			@test isapprox(d[-1,-2],-sin(β)*(1+cos(β))/2,atol=1e-14,rtol=1e-8)

			@test isapprox(d[-2,2],(1-cos(β))^2/4,atol=1e-14,rtol=1e-8)
			@test isapprox(d[-2,-2],(1+cos(β))^2/4,atol=1e-14,rtol=1e-8)
		end

	    n = 100
		j = 2
		checkdjvalues(test,j,n)
	end
	@testset "shift angles" begin
		@testset "β + 4π*n" begin
			for j in 0:1//2:2, β in -4π:π/4:4π
		    	d1 = WignerdMatrix(j, β)
		    	d2 = WignerdMatrix(j, β + 4π)
		    	d3 = WignerdMatrix(j, β - 4π)
		    	@test collect(d1) ≈ collect(d2)
		    	@test collect(d1) ≈ collect(d3)
		    end
		end
		@testset "β + (2n+1)*2π" begin
		    for j in 0:1//2:2, β in -4π:π/2:4π
		    	d1 = WignerdMatrix(j, β)
		    	d2 = WignerdMatrix(j, β + 2π)
		    	d3 = WignerdMatrix(j, β - 2π)
		    	d4 = WignerdMatrix(j, β + 6π)
		    	d5 = WignerdMatrix(j, β - 6π)
		    	@test collect(d1) ≈ (-1)^2j .* collect(d2)
		    	@test collect(d1) ≈ (-1)^2j .* collect(d3)
		    	@test collect(d1) ≈ (-1)^2j .* collect(d4)
		    	@test collect(d1) ≈ (-1)^2j .* collect(d5)
		    end 
		end
		@testset "negative β" begin
		    for j in 0:1//2:2, β in -4π:π/4:4π
		    	d1 = WignerdMatrix(j, β)
		    	d2 = WignerdMatrix(j, -β)
		    	for m in axes(d1,1), n in axes(d1,2)
		    		@test isapprox(d1[m,n],d2[n,m],atol=1e14,rtol=1e-8)
		    	end
		    end
		end
		@testset "π - β" begin
		    for j in 0:1//2:2, β in -4π:π/4:4π
		    	d1 = WignerdMatrix(j, β)
		    	d2 = WignerdMatrix(j, π - β)
		    	for m in axes(d1,1), n in axes(d1,2)
		    		@test isapprox(d2[m,n],(-1)^(j-n) * d1[-m,n],atol=1e-14,rtol=1e-8)
		    	end
		    end
		end
		@testset "β + (2n+1)π" begin
		    for j in 0:1//2:2, β in -4π:π/4:-4π, k = -1:1
		    	d1 = WignerdMatrix(j, β)
		    	d2 = WignerdMatrix(j, β + (2k+1)π)
		    	d3 = WignerdMatrix(j, β - (2k+1)π)
		    	for m in axes(d1,1), n in axes(d1,2)
		    		@test isapprox(d2[m,n],(-1)^((2k+1)j-n) * d1[m,-n],atol=1e-14,rtol=1e-8)
		    		@test isapprox(d3[m,n],(-1)^(-(2k+1)j-n) * d1[m,-n],atol=1e-14,rtol=1e-8)
		    	end
		    end
		end
	end
	@testset "similar" begin
	    d = WignerdMatrix{Float32}(1, π/3)
	    d′ = similar(d)
	    @test d′ isa WignerdMatrix{Float32}
	    @test sphericaldegree(d′) == sphericaldegree(d)
	    @test d′.beta == d.beta
	    @test typeof(d′.dj) == typeof(d.dj)
	    @test size(d′.dj) == size(d.dj)

	    d′ = similar(d, Float64)
	    @test d′ isa WignerdMatrix{Float64}
	    @test sphericaldegree(d′) == sphericaldegree(d)
	    @test d′.beta == d.beta
	    @test eltype(d′) == Float64
	    @test size(d′.dj) == size(d.dj)
	end
	@testset "LinearAlgebra" begin
		@testset "product" begin
			β = π/3
			for j in [1/2,1,3/2,2]
				d = WignerdMatrix(j,β)
				@test collect(d*d) ≈ collect(d)*collect(d)
			end
		end
		@testset "symmetry" begin
			for j = [1/2, 1]
				d = WignerdMatrix(j ,π/3)
			    @test !issymmetric(d)
				@test !ishermitian(d)
			end
		end
		@testset "det" begin
			for j = [1/2, 1]
				d = WignerdMatrix(j ,π/3)
		    	@test det(d) == 1
		    end
		end
		@testset "trace" begin
			for β in LinRange(-4π, 4π, 100)
			    d = WignerdMatrix(1/2,β)
			    @test isapprox(tr(d), 2cos(β/2), atol=1e-14, rtol=1e-8)
			    @test isapprox(tr(d),sum(d[m,m] for m=-1//2:1//2),atol=1e-14,rtol=1e-8)

			    d = WignerdMatrix(1,β)
			    @test isapprox(tr(d) , 1 + 2cos(β), atol=1e-14, rtol=1e-8)
			    @test isapprox(tr(d),sum(d[m,m] for m=-1:1),atol=1e-14,rtol=1e-8)
			    @test isapprox(tr(d),tr(collect(d)),atol=1e-14,rtol=1e-8)

			    d = WignerdMatrix(3/2,β)
			    @test isapprox(tr(d) , 2(cos(β/2) + cos(3β/2)), atol=1e-14, rtol=1e-8)
			    @test isapprox(tr(d),sum(d[m,m] for m=-3//2:3//2),atol=1e-14,rtol=1e-8)
			    @test isapprox(tr(d),tr(collect(d)),atol=1e-14,rtol=1e-8)

			    d = WignerdMatrix(2,β)
			    @test isapprox(tr(d),1 + 2cos(β) + 2cos(2β),atol=1e-14,rtol=1e-8)
			    @test isapprox(tr(d),sum(d[m,m] for m=-2:2),atol=1e-14,rtol=1e-8)
			    @test isapprox(tr(d), tr(collect(d)),atol=1e-14,rtol=1e-8)
			end

			for n = 0:10
			    d1 = WignerdMatrix(n/2, ZeroRadians())
			    d2 = WignerdMatrix(n/2, TwoNPiRadians(0))
			    d3 = WignerdMatrix(n/2, TwoNPiRadians(2))
			    @test tr(d1) == tr(d2) == tr(d3) == Float64(n + 1)
				
				d1 = WignerdMatrix(n/2, π)
				d2 = WignerdMatrix(n/2, TwoNPlusOnePiRadians(0))
				rexp = isinteger(n/2) ? (-1)^Int(n/2) : 0
			    @test tr(d1) == tr(d2) == rexp

			    d = WignerdMatrix(n/2, TwoNPiRadians(1))
			    rexp = isinteger(n/2) ? n+1 : -(n+1)
			    @test tr(d) == rexp
			end
		end
		@testset "inv" begin
			rotangles = -4π:π/4:4π
			js = 0:1//2:4
			for β = rotangles, j in js
				d = WignerdMatrix(j,β)
				dinv = inv(d)
				@test dinv isa WignerdMatrix
				@test collect(dinv) ≈ inv(collect(d))
				@test collect(dinv * d) ≈ Diagonal(ones(size(d,1)))
			end

			for j in js
				d = WignerdMatrix(j, ZeroRadians())
				@test inv(d) == d
				@test inv(d).beta === ZeroRadians()

				d = WignerdMatrix(j, TwoNPlusOnePiRadians(0))
				dinv = WignerdMatrix(j, -TwoNPlusOnePiRadians(0))
				@test inv(d) == dinv
				@test inv(d).beta === -TwoNPlusOnePiRadians(0)

				d = WignerdMatrix(j, TwoNPiRadians(1))
				dinv = WignerdMatrix(j, -TwoNPiRadians(1))
				@test inv(d) == dinv
				@test inv(d).beta === -TwoNPiRadians(1)

				d = WignerdMatrix(j, TwoNPlusOnePiRadians(1))
				dinv = WignerdMatrix(j, -TwoNPlusOnePiRadians(1))
				@test inv(d) == dinv
				@test inv(d).beta === -TwoNPlusOnePiRadians(1)
			end
		end
		@testset "diag" begin
			function testdiag(d)
				j = sphericaldegree(d)
				@test begin 
					res = diag(d) == [d[m,m] for m=-j:j] == diag(collect(d))
					if !res
						@show j, d.beta
					end
					res
				end
			end
			for j = 0:1//2:10
		    	for β in LinRange(-4π, 4π, 100)
			    	d = WignerdMatrix(j, β)
			    	testdiag(d)
			    end

			    for T in [ZeroRadians, Piby2Radians]
			    	d = WignerdMatrix(j, T())
			    	testdiag(d)
			    end

			    d = WignerdMatrix(j, π)
			    testdiag(d)

			    for T in [TwoNPiRadians, TwoNPlusOnePiRadians], n=-3:3
			    	d = WignerdMatrix(j, T(n))
			    	testdiag(d)
			    end
		    end
		end
	end
end
@testset "WignerDMatrix" begin
    @testset "Constructor" begin
    	β = π/3
		α,γ = π/6, π/2
		j = 2
		D = WignerDMatrix(j,α,β,γ)

        @test D.alpha == α
        @test D.gamma == γ
        @test D.dj.beta == β
        @test WignerDMatrices.eulerangles(D) == (α,β,γ)
        @test D.dj isa WignerdMatrix
        @test D.dj == WignerdMatrix(j,β)
        @test WignerDMatrices.sphericaldegree(D) == j

        D2 = WignerDMatrix(j,α,β,γ)
        @test D == D2

        D3 = WignerDMatrix(j,0,β,0)
        d = WignerdMatrix(j,β)
        @test D3 == d
        @test d == D3

        D4 = WignerDMatrix(α, d, γ)
        @test D4.dj === d

        @testset "special angles" begin
        	β = π/3
			α,γ = π/6, π/2
			j = 2
            D = WignerDMatrix(j,TwoNPlusOnePiRadians(0),β,γ)
            @test D[1,1] == -D.dj[1,1] * cis(-γ)
            D = WignerDMatrix(j,TwoNPlusOnePiRadians(0),β,ZeroRadians())
            @test D[1,1] == -D.dj[1,1]
            D = WignerDMatrix(j,α,β,TwoNPlusOnePiRadians(0))
            @test D[1,1] == -D.dj[1,1] * cis(-α)
        end
    end
    @testset "indexing" begin
    	β = π/3
		α,γ = π/6, π/2
		j = 2
		D = WignerDMatrix(j,α,β,γ)
        @test D[1,1] ≈ D.dj[1,1]*cis(-(D.alpha + D.gamma))

        for n in axes(D,2), m in axes(D,1)
        	Dmn = WignerDMatrices.wignerDmatrixelement(j, m, n, (α,β,γ))
        	@test isapprox(D[m,n], Dmn, atol=1e-13, rtol=1e-8)
        end
    end
    @testset "similar" begin
	    D = WignerDMatrix{ComplexF64}(1, 0, π/3, π/3)
	    D′ = similar(D)
	    @test D′ isa WignerDMatrix{ComplexF64}
	    @test sphericaldegree(D′.dj) == sphericaldegree(D.dj)
	    @test D′.dj.beta == D.dj.beta
	    @test typeof(D′.dj) == typeof(D.dj)
	    @test size(D′.dj) == size(D.dj)

	    D′ = similar(D, ComplexF32)
	    @test D′ isa WignerDMatrix{ComplexF32}
	    @test sphericaldegree(D′.dj) == sphericaldegree(D.dj)
	    @test typeof(D′.dj) == typeof(D.dj)
	    @test size(D′.dj) == size(D.dj)
	end
	@testset "LinearAlgebra" begin
		rotangles = -4π:π/3:4π
		@testset "product" begin
			β = π/3
			α,γ = π/6, π/2
			for j in [1/2,1,3/2,2]
				D = WignerDMatrix(j,α,β,γ)
				@test collect(D*D) == collect(D)*collect(D)
				@test collect(D*D.dj) ≈ collect(D)*collect(D.dj)
				@test collect(D.dj*D) ≈ collect(D.dj)*collect(D)
			end
		end
		@testset "symmetry" begin
			for j = [1/2, 1]
				D = WignerDMatrix(j, π/2, π/3, π/4)
			    @test !issymmetric(D)
				@test !ishermitian(D)
			end
		end
		@testset "det" begin
			for j in 0:1//2:4, α in rotangles, 
				β in rotangles, γ in rotangles

				D = WignerDMatrix(j, α, β, γ)
		    	@test det(D) == 1
		    end
		end
		@testset "trace" begin
			for j in 0:1//2:10, α in rotangles, 
				β in rotangles, γ in rotangles
			    D = WignerDMatrix(j,α,β,γ)
			    @test begin 
			    	res = isapprox(tr(D),sum(D[m,m] for m=-j:j),atol=1e-12,rtol=1e-8)
			    	if !res
			    		@show j, α, β, γ
			    	end
			    	res
			    end
			    @test isapprox(tr(D),tr(collect(D)),atol=1e-12,rtol=1e-8)
			end

			for j in 0:1//2:10, α in rotangles, γ in rotangles
				D1 = WignerDMatrix(j,α,ZeroRadians(),γ)
				D2 = WignerDMatrix(j,α,0,γ)
				D3 = WignerDMatrix(j,α,TwoNPiRadians(0),γ)
			    @test isapprox(tr(D1),sum(D1[m,m] for m=-j:j),atol=1e-12,rtol=1e-8)
			    @test isapprox(tr(D1),tr(collect(D1)),atol=1e-12,rtol=1e-8)
			    @test tr(D1) == tr(D2) == tr(D3)

			    D1 = WignerDMatrix(j,α,TwoNPlusOnePiRadians(0),γ)
			    D2 = WignerDMatrix(j,α,π,γ)
			    D3 = WignerDMatrix(j,α,TwoNPlusOnePiRadians(0),γ)
			    @test isapprox(tr(D1),sum(D1[m,m] for m=-j:j),atol=1e-12,rtol=1e-8)
			    @test isapprox(tr(D1),tr(collect(D1)),atol=1e-12,rtol=1e-8)
			    @test tr(D1) == tr(D2) == tr(D3)
			end
		end
		@testset "inv" begin
			for j in 0:1//2:4, α in rotangles, 
				β in rotangles, γ in rotangles

				D = WignerDMatrix(j,α,β,γ)
				Dinv = inv(D)
				@test Dinv isa WignerDMatrix
				@test collect(Dinv) ≈ inv(collect(D))
			end
		end
		@testset "diag" begin
		    function testdiag(d)
				j = sphericaldegree(d)
				@test begin 
					res = diag(d) == [d[m,m] for m=-j:j] == diag(collect(d))
					if !res
						@show j, d.beta
					end
					res
				end
			end
			for j = 0:1//2:10
		    	for α in rotangles, γ in rotangles

		    		for β in rotangles
				    	D = WignerDMatrix(j, α, β, γ)
				    	testdiag(D)
				    end
				    for T in [ZeroRadians, Piby2Radians]
				    	D = WignerDMatrix(j, α, T(), γ)
				    	testdiag(D)
				    end

				    D = WignerDMatrix(j, α, π, γ)
				    testdiag(D)

				    for T in [TwoNPiRadians, TwoNPlusOnePiRadians], n=-3:3
				    	D = WignerDMatrix(j, α, T(n), γ)
				    	testdiag(D)
				    end
			    end
		    end
		end
	end
end
@testset "subarray" begin
	j = half(1); α, β, γ = 2π*rand(), 2π*rand(), 2π*rand()

	function testviewindexing(d)
		dv = @view d[0.5,0.5]
        @test dv[] == d[0.5,0.5]

        dv = @view d[:,0.5]
        for i in axes(d,1)
			@test dv[i] == d[i,0.5]
		end

		dv = @view d[-0.5:0.5, 0.5]
        for (ind,i) in enumerate(axes(d,1))
			@test dv[ind] == d[i,0.5]
		end

		dv = @view d[:,:]
		for i in eachindex(d)
			@test dv[i] == d[i]
		end
	end
    @testset "WignerdMatrix" begin
        d = WignerdMatrix(j, β)
        testviewindexing(d)
    end
    @testset "WignerDMatrix" begin
        d = WignerDMatrix(j, α, β, γ)
        testviewindexing(d)
    end
end
@testset "show" begin
    io = IOBuffer()

    function testshow(io, d)
    	summary(io,d)
	    show(io,MIME"text/plain"(), d)
	    show(io, d)
    end

    for β in [rand(), ZeroRadians(), Piby2Radians(), 
    	TwoNPiRadians(1), TwoNPlusOnePiRadians(1)]
	    
	    d = WignerdMatrix(1, β)
	    testshow(io, d)

	    D = WignerDMatrix(0, d, 0)
    	testshow(io, D)

	    d = WignerdMatrix(0.5, β)
	    testshow(io, d)

	    D = WignerDMatrix(0, d, 0)
    	testshow(io, D)
	end    

    D = WignerDMatrix(1/2, 0, 0, 0)
	testshow(io, D)    

    Dv = @view D[:,:]
    testshow(io, Dv)
end
