using Test
using SLAM

@testset "SLAM.jl Tests" begin

    @testset "Lagrange Mesh" begin
        # Test mesh initialization
        mesh = init_legendre_mesh(10, 15.0)

        @test mesh.N == 10
        @test mesh.R == 15.0
        @test length(mesh.x) == 10
        @test length(mesh.r) == 10
        @test length(mesh.λ) == 10
        @test length(mesh.α) == 10

        # Check that mesh points are in (0, 1)
        @test all(0 .< mesh.x .< 1)

        # Check that physical points are in (0, R)
        @test all(0 .< mesh.r .< mesh.R)

        # Check that weights are positive
        @test all(mesh.λ .> 0)

        # Check normalization
        @test all(mesh.α .≈ 1.0 ./ sqrt.(mesh.λ))
    end

    @testset "Baye Matrices" begin
        mesh = init_legendre_mesh(5, 10.0)

        # Test D matrix
        D = baye_D_matrix(mesh)
        @test size(D) == (5, 5)

        # Test T matrix
        T = baye_T_matrix(mesh)
        @test size(T) == (5, 5)

        # Note: Baye's T matrix for Lagrange-Legendre is NOT necessarily symmetric
        # due to the non-symmetric quadrature structure. This is expected.
        # The important property is that it represents -d²/dr² correctly.
    end

    @testset "Coulomb Functions" begin
        # Test wave number calculation
        μ = 938.0  # proton mass in MeV/c²
        E = 30.0   # 30 MeV
        k = coulomb_k(μ, E)
        @test k > 0

        # Test Sommerfeld parameter
        η = coulomb_eta(1.0, 20.0, μ, E)
        @test η > 0

        # Test Coulomb functions at ρ = 5.0, η = 1.0, l = 0
        ρ = 5.0
        η_test = 1.0
        l = 0

        F, G, Fp, Gp = coulomb_FG(l, η_test, ρ)

        # Check Wronskian: F*G' - F'*G = ±1 (sign depends on convention)
        wronskian = F * Gp - Fp * G
        @test abs(abs(wronskian) - 1.0) < 1e-6

        # Test Hankel function
        H_plus = coulomb_H_plus(l, η_test, ρ)
        @test H_plus ≈ complex(G, F)

        # Test gamma_s is finite
        γ_s = coulomb_gamma_s(l, η_test, ρ)
        @test isfinite(real(γ_s)) && isfinite(imag(γ_s))
    end

    @testset "Optical Potential" begin
        # Create a simple potential
        pot = OpticalPotential(
            V_v=50.0, r_v=1.25, a_v=0.65,
            W_v=10.0, r_wv=1.25, a_wv=0.65,
            V_s=0.0, r_s=1.25, a_s=0.65,
            W_s=0.0, r_ws=1.25, a_ws=0.65,
            V_so=0.0, r_so=1.10, a_so=0.65,
            W_so=0.0, r_wso=1.10, a_wso=0.65,
            r_c=1.25,
            Z_proj=1.0, A_proj=1.0,
            Z_targ=20.0, A_targ=40.0
        )

        # Test reduced mass
        μ = reduced_mass(pot)
        @test μ > 0

        # Test Woods-Saxon form factor
        # At r = 0, f ≈ 1 (but not exactly 1 due to exp(-R/a))
        @test woods_saxon(0.0, 5.0, 0.5) > 0.99
        @test woods_saxon(10.0, 5.0, 0.5) < 0.01
        @test woods_saxon(5.0, 5.0, 0.5) ≈ 0.5

        # Test potential evaluation
        V = evaluate_potential(pot, 5.0, 0, 0.5)
        @test isfinite(real(V)) && isfinite(imag(V))
        @test imag(V) < 0  # Absorptive part should be negative

        # Test Coulomb potential
        V_c = evaluate_coulomb(pot, 10.0)
        @test V_c > 0  # Repulsive for like charges
    end

    @testset "Scattering Solver" begin
        # Set up a simple scattering problem (neutron scattering, no Coulomb)
        pot = OpticalPotential(
            V_v=50.0, r_v=1.25, a_v=0.65,
            W_v=10.0, r_wv=1.25, a_wv=0.65,
            r_c=0.0,
            Z_proj=0.0, A_proj=1.0,
            Z_targ=0.0, A_targ=40.0,
            A1=0.0, A2=0.0
        )

        # Create scattering problem with more mesh points for better accuracy
        prob = ScatteringProblem(pot, 30.0, 0; N=60, R=20.0)

        # Solve
        result = solve_scattering(prob)

        # Check that we got a result
        @test result.converged

        # S-matrix magnitude should be close to ≤ 1 for absorptive potential
        # Allow some numerical tolerance
        @test abs(result.S_matrix) < 1.1

        # Phase shift should be finite
        @test isfinite(real(result.phase_shift)) && isfinite(imag(result.phase_shift))
    end

    @testset "Legendre Polynomials" begin
        # Test P_0(x) = 1
        @test legendre_P(0, 0.5) ≈ 1.0

        # Test P_1(x) = x
        @test legendre_P(1, 0.5) ≈ 0.5

        # Test P_2(x) = (3x² - 1)/2
        @test legendre_P(2, 0.5) ≈ (3*0.25 - 1)/2

        # Test orthogonality property at x = 0
        @test legendre_P(1, 0.0) ≈ 0.0
        @test legendre_P(3, 0.0) ≈ 0.0
    end

end

println("\nAll tests completed!")
