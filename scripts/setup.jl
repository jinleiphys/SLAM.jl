#!/usr/bin/env julia
"""
    SLAM.jl Setup Script

Installs all required dependencies for the SLAM.jl GUI.

Usage:
    julia setup.jl
"""

using Pkg

println("="^60)
println("    SLAM.jl - Setup Script")
println("="^60)
println()

# Activate project
println("Activating project...")
Pkg.activate(@__DIR__)

# Install dependencies
println("\nInstalling dependencies...")
println("-"^60)

deps = [
    "Blink",
    "JSON3",
    "FastGaussQuadrature",
    "SpecialFunctions",
    "LinearAlgebra",
    "Printf"
]

for dep in deps
    print("  Installing $dep... ")
    try
        Pkg.add(dep)
        println("OK")
    catch e
        println("FAILED: $e")
    end
end

# Instantiate to resolve versions
println("\nResolving package versions...")
Pkg.instantiate()

# Precompile
println("\nPrecompiling packages...")
Pkg.precompile()

# Build Blink (downloads Electron)
println("\nBuilding Blink (downloading Electron)...")
try
    Pkg.build("Blink")
    println("  Blink build complete!")
catch e
    println("  Note: Blink build issue (may still work): $e")
end

println()
println("="^60)
println("    Setup Complete!")
println("="^60)
println()
println("To start the GUI, run:")
println("    julia run_gui.jl")
println()
println("Or from Julia REPL:")
println("    include(\"gui.jl\")")
println("    main()")
println()
