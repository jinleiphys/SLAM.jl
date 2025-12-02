#!/usr/bin/env julia
"""
    SLAM.jl GUI Launcher

One-click launcher for the SLAM.jl graphical interface.

Usage:
    julia run_gui.jl
"""

# Activate project (root directory)
using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

# Install dependencies if needed
println("Checking dependencies...")

# List of required packages for GUI
gui_packages = ["Blink", "JSON3"]

# Check and install each package
for pkg in gui_packages
    if !haskey(Pkg.project().dependencies, pkg)
        println("Installing $pkg...")
        Pkg.add(pkg)
    end
end

# Instantiate to ensure everything is ready
Pkg.instantiate()

println("Dependencies OK")

# Run GUI
include(joinpath(@__DIR__, "gui.jl"))
main()
