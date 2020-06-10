using WignerDMatrices
using HalfIntegers
using LinearAlgebra
using Test

import WignerDMatrices: cis_special, BetaPiby2, BetaZero, BetaPi

@testset "special points" begin
	@testset "float" begin
		@testset "BetaPiby2" begin
		    @test float(BetaPiby2()) == π/2
	    	@test AbstractFloat(BetaPiby2()) == π/2
	    	@test Float64(BetaPiby2()) == π/2
		end
		@testset "BetaZero" begin
			@test float(BetaZero()) == 0
	    	@test AbstractFloat(BetaZero()) == 0
	    	@test Float64(BetaZero()) == 0
		end
		@testset "BetaPi" begin
			@test float(BetaPi()) == float(π)
	    	@test AbstractFloat(BetaPi()) == float(π)
	    	@test Float64(BetaPi()) == float(π)
		end
	end
	@testset "promote" begin
		for T in [Int,Float32,Float64]
		    @test promote_rule(BetaPiby2,T) == Float64
		    @test promote_rule(BetaZero,T) == Float64
		    @test promote_rule(BetaPi,T) == Float64
		end
		for T in [BigFloat]
			@test promote_rule(BetaPiby2,T) == BigFloat
			@test promote_rule(BetaZero,T) == BigFloat
		    @test promote_rule(BetaPi,T) == BigFloat
		end
	end
	@testset "zero and one" begin
		@testset "BetaPiby2" begin
		    @test zero(BetaPiby2()) == zero(Float64)
	    	@test one(BetaPiby2()) == one(Float64)
		end
		@testset "BetaZero" begin
    		@test zero(BetaZero()) == zero(Float64)
    		@test one(BetaZero()) == one(Float64)
		end
		@testset "BetaPi" begin
	    	@test zero(BetaPi()) == zero(Float64)
	    	@test one(BetaPi()) == one(Float64)
		end
	end
	@testset "trigonometric functions" begin
	    @testset "BetaPiby2" begin
			@test one(BetaPiby2()) == 1
		    for α = -10:1//2:10
		    	@test cis_special(HalfInt(α),BetaPiby2()) ≈ cis(α*π/2)
		    end
	    	@test cos(BetaPiby2()) == 0
	    	@test sin(BetaPiby2()) == 1
		end
		@testset "North Pole" begin
			for α = -10:1//2:10
		    	@test cis_special(HalfInt(α),BetaZero()) ≈ cis(0)
		    end
		    @test cos(BetaZero()) == 1
	    	@test sin(BetaZero()) == 0
		end
		@testset "South Pole" begin
			for α = -10:1//2:10
		    	@test cis_special(HalfInt(α),BetaPi()) ≈ cis(α*π)
		    end
		    @test cos(BetaPi()) == -1
	    	@test sin(BetaPi()) == 0
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
	    	@testset "BetaPiby2" begin
		        for m in -j:j, n in -j:j
		        	dj_m_n = WignerDMatrices.djmatrix_terms(π/2,λ,v,m,n,j)
		        	dj_m_n2 = WignerDMatrices.djmatrix_terms(BetaPiby2(),λ,v,m,n,j)

		        	testapprox(m,n,dj_m_n,dj_m_n2)
		        end
		    end
		    @testset "BetaZero" begin
		        for m in -j:j, n in -j:j
		        	dj_m_n = WignerDMatrices.djmatrix_terms(0,λ,v,m,n,j)
		        	dj_m_n2 = WignerDMatrices.djmatrix_terms(BetaZero(),λ,v,m,n,j)

		        	testapprox(m,n,dj_m_n,dj_m_n2)
		        end
		    end
		    @testset "BetaPi" begin
		        for m in -j:j, n in -j:j
		        	dj_m_n = WignerDMatrices.djmatrix_terms(π,λ,v,m,n,j)
		        	dj_m_n2 = WignerDMatrices.djmatrix_terms(BetaPi(),λ,v,m,n,j)

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
			@test axes(d) == (-jh:jh,-jh:jh)
			@test axes(d,1) == axes(d,2) == -jh:jh
			@test size(d) == (Integer(2j+1),Integer(2j+1))
			@test size(d,1) == size(d,2) == Integer(2j+1)

			@test d == d

			d = WignerdMatrix(j, zeros(WignerDMatrices.filledelements(j)+1))
			@test WignerDMatrices.sphericaldegree(d) == j
			@test axes(d) == (-jh:jh,-jh:jh)

			d = WignerdMatrix(j, zeros(WignerDMatrices.filledelements(j)))
			@test WignerDMatrices.sphericaldegree(d) == j
			@test axes(d) == (-jh:jh,-jh:jh)

			@test_throws ArgumentError WignerdMatrix(j, zeros(0))

			for i in eachindex(d.dj)
				@test d[i] == d.dj[i]
			end
		end
		@testset "half-integer j" begin
			j = 3//2
			jh = HalfInt(j)
			d = WignerdMatrix(j,β)
			@test WignerDMatrices.sphericaldegree(d) == j
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

		@testset "BetaPiby2" begin
			β = BetaPiby2()
			d = WignerdMatrix(j,β)
			test(d,β) 
		end
		@testset "BetaZero" begin
			β = BetaZero()
			d = WignerdMatrix(j,β)
			test(d,β) 
		end
		@testset "BetaPi" begin
			β = BetaPi()
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
	    @test typeof(d′.dj) == typeof(d.dj)
	    @test size(d′.dj) == size(d.dj)

	    d′ = similar(d, Float64)
	    @test d′ isa WignerdMatrix{Float64}
	    @test d′.j == d.j
	    @test typeof(d′.dj) == typeof(d.dj)
	    @test size(d′.dj) == size(d.dj)
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
    end
    @testset "update angles" begin
        β = π/3
    	α,γ = π/6, π/2
    	j = 2
        D = WignerDMatrix(j,α,β,γ)
        WignerdMatrix!(D,π/2)
        @test D.dj == WignerdMatrix(j,π/2)

        α_new, γ_new = π/12, π/12
        D.α = α_new
        @test D.α == α_new 
        D.γ = γ_new
        @test D.γ == γ_new 
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
	    @test typeof(D′.dj) == typeof(D.dj)
	    @test size(D′.dj) == size(D.dj)

	    D′ = similar(D, ComplexF32)
	    @test D′ isa WignerDMatrix{ComplexF32}
	    @test D′.dj.j == D.dj.j
	    @test typeof(D′.dj) == typeof(D.dj)
	    @test size(D′.dj) == size(D.dj)
	end
end
@testset "show" begin
    io = IOBuffer()
    d = WignerdMatrix(1, π/2)
    summary(io,d)
    @test String(take!(io)) == "Wigner d-matrix with j = $(d.j)"
    D = WignerDMatrix(0, d, 0)
    summary(io,D)
    strexp = "Wigner D-matrix with j = $(d.j)"*
    ", with α = $(D.α) and γ = $(D.γ)"
    @test String(take!(io)) == strexp

    show(io,MIME"text/plain"(), d)
    show(io, d)
end
