# Debug: Compare mesh points
using SLAM
using Printf

function gauleg_fortran(n)
    # Exact Fortran gauleg algorithm
    eps = 3e-14
    x = zeros(n)
    w = zeros(n)

    m = div(n + 1, 2)

    for i in 1:m
        z = cos(π * (i - 0.25) / (n + 0.5))
        local pp
        while true
            p1 = 1.0
            p2 = 0.0
            for j in 1:n
                p3 = p2
                p2 = p1
                p1 = ((2*j - 1) * z * p2 - (j - 1) * p3) / j
            end
            pp = n * (z * p1 - p2) / (z^2 - 1)
            z1 = z
            z = z1 - p1 / pp
            if abs(z - z1) < eps
                break
            end
        end
        x[i] = -z
        x[n + 1 - i] = z
        w[i] = 2.0 / ((1 - z^2) * pp^2)
        w[n + 1 - i] = w[i]
    end

    return x, w
end

function main()
    N = 60
    R = 20.0

    println("=" ^ 70)
    println("Mesh point comparison (N=$N, R=$R)")
    println("=" ^ 70)
    println()

    # Get SLAM mesh
    mesh = init_legendre_mesh(N, R)

    # Get Fortran-style mesh
    t_fort, w_fort = gauleg_fortran(N)
    x_fort = (t_fort .+ 1) ./ 2  # Map to (0,1)
    w_fort = w_fort ./ 2  # Scale weights
    r_fort = R .* x_fort

    println("First 5 mesh points:")
    println("  j      SLAM x         Fortran x       diff(x)           SLAM r         Fortran r")
    println("-" ^ 90)
    for j in 1:5
        diff_x = mesh.x[j] - x_fort[j]
        @printf("  %d   %.12f   %.12f   %+.2e   %.10f   %.10f\n",
                j, mesh.x[j], x_fort[j], diff_x, mesh.r[j], r_fort[j])
    end

    println()
    println("Last 5 mesh points:")
    for j in (N-4):N
        diff_x = mesh.x[j] - x_fort[j]
        @printf("  %d   %.12f   %.12f   %+.2e   %.10f   %.10f\n",
                j, mesh.x[j], x_fort[j], diff_x, mesh.r[j], r_fort[j])
    end

    println()
    println("Weight comparison:")
    println("  j      SLAM λ         Fortran λ       diff")
    println("-" ^ 60)
    for j in 1:5
        diff_w = mesh.λ[j] - w_fort[j]
        @printf("  %d   %.12f   %.12f   %+.2e\n",
                j, mesh.λ[j], w_fort[j], diff_w)
    end

    max_diff_x = maximum(abs.(mesh.x .- x_fort))
    max_diff_w = maximum(abs.(mesh.λ .- w_fort))
    println()
    println("Maximum differences:")
    println("  max|x - x_fort| = $max_diff_x")
    println("  max|λ - λ_fort| = $max_diff_w")
end

main()
