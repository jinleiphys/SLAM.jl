#!/bin/bash
# SLAM.jl Setup Script (Shell Script)
#
# This script will:
#   1. Check if Julia is installed
#   2. If not, install Julia via Homebrew (macOS) or juliaup
#   3. Install all required Julia packages
#
# Usage:
#   ./setup.sh
#
# Or make executable first:
#   chmod +x setup.sh
#   ./setup.sh

set -e

echo "============================================================"
echo "    SLAM.jl - Setup Script"
echo "============================================================"
echo ""

cd "$(dirname "$0")"

# Function to check if a command exists
command_exists() {
    command -v "$1" >/dev/null 2>&1
}

# Check if Julia is installed
if command_exists julia; then
    JULIA_VERSION=$(julia --version 2>&1)
    echo "✓ Julia is installed: $JULIA_VERSION"
else
    echo "Julia is not installed. Installing..."
    echo ""

    # Detect OS
    OS="$(uname -s)"

    case "$OS" in
        Darwin)
            # macOS
            echo "Detected macOS"

            if command_exists brew; then
                echo "Installing Julia via Homebrew..."
                brew install julia
            elif command_exists curl; then
                echo "Installing Julia via juliaup..."
                curl -fsSL https://install.julialang.org | sh -s -- -y
                # Add to PATH for current session
                export PATH="$HOME/.juliaup/bin:$PATH"
            else
                echo "ERROR: Neither Homebrew nor curl found."
                echo "Please install Julia manually from https://julialang.org/downloads/"
                exit 1
            fi
            ;;

        Linux)
            # Linux
            echo "Detected Linux"

            if command_exists curl; then
                echo "Installing Julia via juliaup..."
                curl -fsSL https://install.julialang.org | sh -s -- -y
                # Add to PATH for current session
                export PATH="$HOME/.juliaup/bin:$PATH"
            elif command_exists wget; then
                echo "Installing Julia via juliaup (wget)..."
                wget -qO- https://install.julialang.org | sh -s -- -y
                export PATH="$HOME/.juliaup/bin:$PATH"
            else
                echo "ERROR: Neither curl nor wget found."
                echo "Please install Julia manually from https://julialang.org/downloads/"
                exit 1
            fi
            ;;

        MINGW*|MSYS*|CYGWIN*)
            # Windows (Git Bash, MSYS, Cygwin)
            echo "Detected Windows"
            echo "Please install Julia manually from https://julialang.org/downloads/"
            echo "Or run in PowerShell: winget install julia -s msstore"
            exit 1
            ;;

        *)
            echo "Unknown operating system: $OS"
            echo "Please install Julia manually from https://julialang.org/downloads/"
            exit 1
            ;;
    esac

    # Verify installation
    if command_exists julia; then
        JULIA_VERSION=$(julia --version 2>&1)
        echo ""
        echo "✓ Julia installed successfully: $JULIA_VERSION"
    else
        # Try to find julia in common locations
        if [ -f "$HOME/.juliaup/bin/julia" ]; then
            export PATH="$HOME/.juliaup/bin:$PATH"
            echo "✓ Julia installed at ~/.juliaup/bin/julia"
        else
            echo "ERROR: Julia installation failed."
            echo "Please install Julia manually from https://julialang.org/downloads/"
            exit 1
        fi
    fi
fi

echo ""
echo "------------------------------------------------------------"
echo "Installing Julia packages..."
echo "------------------------------------------------------------"
echo ""

# Run Julia setup script
julia --project=. scripts/setup.jl

echo ""
echo "============================================================"
echo "    Setup Complete!"
echo "============================================================"
echo ""
echo "To start the GUI, run:"
echo "    ./run_gui.sh"
echo ""
echo "Or:"
echo "    julia --project=. run_gui.jl"
echo ""
