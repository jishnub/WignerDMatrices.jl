using WignerDMatrices
using HalfIntegers
using LinearAlgebra
using Test

import WignerDMatrices: cis_special, Piby2Radians, ZeroRadians, PiRadians

@testset "special points" begin
	@testset "float" begin
		@testset "Piby2Radians" begin
		    @test float(Piby2Radians()) == π/2
	    	@test AbstractFloat(Piby2Radians()) == π/2
	    	@test Float64(Piby2Radians()) == π/2
		end
		@testset "ZeroRadians" begin
			@test float(ZeroRadians()) == 0
	    	@test AbstractFloat(ZeroRadians()) == 0
	    	@test Float64(ZeroRadians()) == 0
		end
		@testset "PiRadians" begin
			@test float(PiRadians()) == float(π)
	    	@test AbstractFloat(PiRadians()) == float(π)
	    	@test Float64(PiRadians()) == float(π)
		end
	end
	@testset "promote" begin
		for T in [Int,Float32,Float64]
		    @test promote_rule(Piby2Radians,T) == Float64
		    @test promote_rule(ZeroRadians,T) == Float64
		    @test promote_rule(PiRadians,T) == Float64
		end
		for T in [BigFloat]
			@test promote_rule(Piby2Radians,T) == BigFloat
			@test promote_rule(ZeroRadians,T) == BigFloat
		    @test promote_rule(PiRadians,T) == BigFloat
		end
	end
	@testset "zero and one" begin
		@testset "Piby2Radians" begin
		    @test zero(Piby2Radians()) == zero(Float64)
	    	@test one(Piby2Radians()) == one(Float64)
		end
		@testset "ZeroRadians" begin
    		@test zero(ZeroRadians()) == zero(Float64)
    		@test one(ZeroRadians()) == one(Float64)
		end
		@testset "PiRadians" begin
	    	@test zero(PiRadians()) == zero(Float64)
	    	@test one(PiRadians()) == one(Float64)
		end
	end
	@testset "trigonometric functions" begin
	    @testset "Piby2Radians" begin
			@test one(Piby2Radians()) == 1
		    for α = -10:1//2:10
		    	@test cis_special(HalfInt(α),Piby2Radians()) ≈ cis(α*π/2)
		    end
	    	@test cos(Piby2Radians()) == 0
	    	@test sin(Piby2Radians()) == 1
		end
		@testset "North Pole" begin
			for α = -10:1//2:10
		    	@test cis_special(HalfInt(α),ZeroRadians()) ≈ cis(0)
		    end
		    @test cos(ZeroRadians()) == 1
	    	@test sin(ZeroRadians()) == 0
		end
		@testset "South Pole" begin
			for α = -10:1//2:10
		    	@test cis_special(HalfInt(α),PiRadians()) ≈ cis(α*π)
		    end
		    @test cos(PiRadians()) == -1
	    	@test sin(PiRadians()) == 0
		end
	end
	@testset "djmatrix_terms" begin
		function testapprox(m,n,dj_m_n,dj_m_n2)
	    	@test begin 
	    		res = isapprox(dj_m_n,dj_m_n2,atol=1e-14,rtol=sqrt(eps(Float64)))
	    		if !res
	    			@show m n dj_m_n dj_m_n2
	    		end
	    		res
	    	end
	    end
	    function tests(j,λ,v)
	    	@testset "Piby2Radians" begin
		        for m in -j:j, n in -j:j
		        	dj_m_n = WignerDMatrices.djmatrix_terms(π/2,λ,v,m,n,j)
		        	dj_m_n2 = WignerDMatrices.djmatrix_terms(Piby2Radians(),λ,v,m,n,j)

		        	testapprox(m,n,dj_m_n,dj_m_n2)
		        end
		    end
		    @testset "ZeroRadians" begin
		        for m in -j:j, n in -j:j
		        	dj_m_n = WignerDMatrices.djmatrix_terms(0,λ,v,m,n,j)
		        	dj_m_n2 = WignerDMatrices.djmatrix_terms(ZeroRadians(),λ,v,m,n,j)

		        	testapprox(m,n,dj_m_n,dj_m_n2)
		        end
		    end
		    @testset "PiRadians" begin
		        for m in -j:j, n in -j:j
		        	dj_m_n = WignerDMatrices.djmatrix_terms(π,λ,v,m,n,j)
		        	dj_m_n2 = WignerDMatrices.djmatrix_terms(PiRadians(),λ,v,m,n,j)

		        	testapprox(m,n,dj_m_n,dj_m_n2)
		        end
		    end
		end

		@testset "integer j" begin
		    j = 5
			A = zeros(ComplexF64,2j+1,2j+1)
		    v = WignerDMatrices.Jy_eigen!(j,A)
		    λ = -j:j

		    tests(j,λ,v)
		end
		@testset "half-integer j" begin
		    j = 3//2
			A = zeros(ComplexF64,Int(2j+1),Int(2j+1))
		    v = WignerDMatrices.Jy_eigen!(HalfInt(j),A)
		    λ = -j:j

		    tests(j,λ,v)
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
	@testset "promote_type_phase" begin
	    @test WignerDMatrices.promote_type_phase(Float64) == Float64
	    @test WignerDMatrices.promote_type_phase(Float32) == Float32
	end
	@testset "Constructor" begin
		β = π/3
		@testset "integer j" begin
			j = 2
			jh = HalfInt(j)
			d = WignerdMatrix(j,β)
			@test WignerDMatrices.sphericaldegree(d) == j
			@test d.β == β
			@test WignerDMatrices.eulerangles(d) == (zero(β),β,zero(β))
			@test axes(d) == (-jh:jh,-jh:jh)
			@test axes(d,1) == axes(d,2) == -jh:jh
			@test size(d) == (Integer(2j+1),Integer(2j+1))
			@test size(d,1) == size(d,2) == Integer(2j+1)

			@test d == d

			d = WignerdMatrix(j, 0.0, zeros(WignerDMatrices.filledelements(j)+1))
			@test WignerDMatrices.sphericaldegree(d) == j
			@test axes(d) == (-jh:jh,-jh:jh)

			d = WignerdMatrix(j, 0.0, zeros(WignerDMatrices.filledelements(j)))
			@test WignerDMatrices.sphericaldegree(d) == j
			@test axes(d) == (-jh:jh,-jh:jh)

			@test_throws ArgumentError WignerdMatrix(j, 0.0, zeros(0))

			for i in eachindex(d.dj)
				@test d[i] == d.dj[i]
			end
		end
		@testset "half-integer j" begin
			j = 3//2
			jh = HalfInt(j)
			d = WignerdMatrix(j,β)
			@test WignerDMatrices.sphericaldegree(d) == j
			@test d.β == β
			@test WignerDMatrices.eulerangles(d) == (zero(β),β,zero(β))
			@test axes(d) == (-jh:jh,-jh:jh)
			@test axes(d,1) == axes(d,2) == -jh:jh
			@test size(d) == (Integer(2j+1),Integer(2j+1))
			@test size(d,1) == size(d,2) == Integer(2j+1)

			for i in eachindex(d.dj)
				@test d[i] == d.dj[i]
			end

			d[1] = 4
			@test d[1] ≈ 4
		end
	end
	function checkdjvalues(test,j,n)
		for β in LinRange(π/n,π-π/n,2n+1)
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
		@testset "PiRadians" begin
			β = PiRadians()
			d = WignerdMatrix(j,β)
			test(d,β) 
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
	@testset "similar" begin
	    d = WignerdMatrix(Float32, 1, π/3)
	    d′ = similar(d)
	    @test d′ isa WignerdMatrix{Float32}
	    @test d′.j == d.j
	    @test d′.β == d.β
	    @test typeof(d′.dj) == typeof(d.dj)
	    @test size(d′.dj) == size(d.dj)

	    d′ = similar(d, Float64)
	    @test d′ isa WignerdMatrix{Float64}
	    @test d′.j == d.j
	    @test d′.β == d.β
	    @test typeof(d′.dj) == typeof(d.dj)
	    @test size(d′.dj) == size(d.dj)
	end
	@testset "LinearAlgebra" begin
		@testset "product" begin
			β = π/3
			for j in [1/2,1,3/2,2]
				d = WignerdMatrix(j,β)
				v = rand(Int(2j+1))
				m = rand(Int(2j+1),Int(2j+1))
				@test d*v == collect(d)*v
				@test_throws DimensionMismatch v*d
				@test d*m == collect(d)*m
				@test m*d == m*d
				@test v'*d == v'*collect(d)
				@test transpose(v)*d == transpose(v)*collect(d)
				@test d*d == collect(d)*collect(d)
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
			β = π/4
		    d = WignerdMatrix(1/2,β)
		    @test tr(d) ≈ 2cos(β/2)
		    @test tr(d) ≈ sum(d[m,m] for m=-1//2:1//2)

		    d = WignerdMatrix(1/2,ZeroRadians())
		    @test tr(d) == Float64(2)

		    d = WignerdMatrix(1,β)
		    @test tr(d) ≈ 1 + 2cos(β)
		    @test tr(d) ≈ sum(d[m,m] for m=-1:1)

		    d = WignerdMatrix(1,ZeroRadians())
		    @test tr(d) == Float64(3)

		    d = WignerdMatrix(3/2,β)
		    @test tr(d) ≈ 2(cos(β/2) + cos(3β/2))
		    @test tr(d) ≈ sum(d[m,m] for m=-3//2:3//2)

		    d = WignerdMatrix(3/2,ZeroRadians())
		    @test tr(d) == Float64(4)

		    d = WignerdMatrix(2,β)
		    @test tr(d) ≈ 1 + 2cos(β) + 2cos(2β)
		    @test tr(d) ≈ sum(d[m,m] for m=-2:2)

		    d = WignerdMatrix(2,ZeroRadians())
		    @test tr(d) == Float64(5)
		end
		@testset "inv" begin
			for β = π/3:π/3:4π, j in 1//2:1//2:2
				d = WignerdMatrix(j,β)
				dinv = inv(d)
				@test dinv isa WignerdMatrix
				@test dinv * d ≈ I
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

        @test D.α == α
        @test D.γ == γ
        @test D.dj.β == β
        @test WignerDMatrices.eulerangles(D) == (α,β,γ)
        @test D.dj isa WignerdMatrix
        @test D.dj == WignerdMatrix(j,β)
        @test WignerDMatrices.sphericaldegree(D) == j

        @test collect(D) == D
        @test D == collect(D)

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
            D = WignerDMatrix(j,PiRadians(),β,γ)
            @test D[1,1] == -D.dj[1,1] * cis(-γ)
            D = WignerDMatrix(j,PiRadians(),β,ZeroRadians())
            @test D[1,1] == -D.dj[1,1]
            D = WignerDMatrix(j,α,β,PiRadians())
            @test D[1,1] == -D.dj[1,1] * cis(-α)
        end
    end
    @testset "update angles" begin
    	β = π/3
		α,γ = π/6, π/2
		j = 2
		D = WignerDMatrix(j,α,β,γ)
		β_new = π/2
        WignerdMatrix!(D,β_new)
        @test D.dj == WignerdMatrix(j,β_new)
    end
    @testset "indexing" begin
    	β = π/3
		α,γ = π/6, π/2
		j = 2
		D = WignerDMatrix(j,α,β,γ)
        @test D[1,1] == D.dj[1,1]*cis(-(D.α + D.γ))

		D[1,1] = 4*cis(-(D.α + D.γ))
		@test D.dj[1,1] ≈ 4
		@test D[1,1] ≈ 4*cis(-(D.α + D.γ))
    end
    @testset "similar" begin
	    D = WignerDMatrix(ComplexF64, 1, 0, π/3, π/3)
	    D′ = similar(D)
	    @test D′ isa WignerDMatrix{ComplexF64}
	    @test D′.dj.j == D.dj.j
	    @test D′.dj.β == D.dj.β
	    @test typeof(D′.dj) == typeof(D.dj)
	    @test size(D′.dj) == size(D.dj)

	    D′ = similar(D, ComplexF32)
	    @test D′ isa WignerDMatrix{ComplexF32}
	    @test D′.dj.j == D.dj.j
	    @test typeof(D′.dj) == typeof(D.dj)
	    @test size(D′.dj) == size(D.dj)
	end
	@testset "LinearAlgebra" begin
		@testset "product" begin
			β = π/3
			α,γ = π/6, π/2
			for j in [1/2,1,3/2,2]
				D = WignerDMatrix(j,α,β,γ)
				v = rand(Int(2j+1))
				m = rand(Int(2j+1),Int(2j+1))
				@test D*v == collect(D)*v
				@test_throws DimensionMismatch v*D
				@test D*m == collect(D)*m
				@test m*D == m*D
				@test v'*D == v'*collect(D)
				@test transpose(v)*D == transpose(v)*collect(D)
				@test D*D == collect(D)*collect(D)
				@test D*D.dj == collect(D)*collect(D.dj)
				@test D.dj*D == collect(D.dj)*collect(D)
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
			for j = [1/2, 1]
				D = WignerDMatrix(j, π/2, π/3, π/4)
		    	@test det(D) == 1
		    end
		end
		@testset "trace" begin
			β = π/3
			α,γ = π/6, π/2
			for j in 1//2:1//2:2
			    D = WignerDMatrix(j,α,β,γ)
			    @test tr(D) ≈ sum(D[m,m] for m=-j:j)
			end
		end
		@testset "inv" begin
			α,γ = π/6, π/2
			for β = π/3:π/3:4π, j in 1//2:1//2:2
				D = WignerDMatrix(j,α,β,γ)
				Dinv = inv(D)
				@test Dinv isa WignerDMatrix
				@test Dinv * D ≈ I
			end
		end
	end
end
@testset "show" begin
    io = IOBuffer()

    d = WignerdMatrix(1, π/2)
    summary(io,d)
    strexp = "Wigner d-matrix with j = $(d.j) and β = $(d.β)"
    @test String(take!(io)) == strexp

    D = WignerDMatrix(0, d, 0)
    summary(io,D)
    α,β,γ = WignerDMatrices.eulerangles(D)
    strexp = "Wigner D-matrix with j = $(d.j)"*
    ", with α = $α, β = $β and γ = $γ"
    @test String(take!(io)) == strexp

    show(io,MIME"text/plain"(), d)
    show(io, d)
end
