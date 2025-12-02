# Debug: Compare T matrix between SLAM.jl and Baye formula
using SLAM

function main()
    N = 5
    R = 20.0

    mesh = init_legendre_mesh(N, R)
    T = baye_T_matrix(mesh)

    println("=" ^ 60)
    println("T matrix comparison (N=$N, R=$R)")
    println("=" ^ 60)
    println()

    println("Mesh points x (in unit interval):")
    for j in 1:N
        println("  x[$j] = $(mesh.x[j])")
    end
    println()

    println("SLAM.jl T matrix (scaled by R² to get unit interval T):")
    T_unit = T * R^2  # Convert back to unit interval
    for i in 1:N
        for j in 1:N
            print("  $(round(T_unit[i,j], digits=6))")
        end
        println()
    end
    println()

    # Now compute using the exact Baye formula directly
    println("Direct Baye formula (on unit interval):")
    x = mesh.x
    for i in 1:N
        xi = x[i]
        for j in 1:N
            xj = x[j]
            if i == j
                num = N * (N + 1) * xi * (1 - xi) - 3 * xi + 1
                den = 3 * xi^2 * (1 - xi)^2
                T_ij = num / den
            else
                sign_ij = (-1)^(i - j)
                num = xi + xj - 2 * xi^2
                den = xj * (xj - xi)^2
                sqrt_num = xj * (1 - xj)
                sqrt_den = xi * (1 - xi)^3
                T_ij = sign_ij * (num / den) * sqrt(sqrt_num / sqrt_den)
            end
            print("  $(round(T_ij, digits=6))")
        end
        println()
    end
    println()

    # Check if they match
    println("Difference (SLAM - Direct):")
    max_diff = 0.0
    for i in 1:N
        for j in 1:N
            xi = x[i]
            xj = x[j]
            if i == j
                num = N * (N + 1) * xi * (1 - xi) - 3 * xi + 1
                den = 3 * xi^2 * (1 - xi)^2
                T_ij = num / den
            else
                sign_ij = (-1)^(i - j)
                num = xi + xj - 2 * xi^2
                den = xj * (xj - xi)^2
                sqrt_num = xj * (1 - xj)
                sqrt_den = xi * (1 - xi)^3
                T_ij = sign_ij * (num / den) * sqrt(sqrt_num / sqrt_den)
            end
            diff = T_unit[i,j] - T_ij
            max_diff = max(max_diff, abs(diff))
            print("  $(round(diff, digits=10))")
        end
        println()
    end
    println("\nMax difference: $max_diff")
end

main()
