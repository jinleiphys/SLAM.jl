#!/bin/bash
# SLAM.jl GUI Launcher (Shell Script)
#
# This script will:
#   1. Check if Julia is installed
#   2. If not, run setup.sh first
#   3. Launch the GUI
#
# Usage:
#   ./run_gui.sh

cd "$(dirname "$0")"

# Function to check if a command exists
command_exists() {
    command -v "$1" >/dev/null 2>&1
}

# Check if Julia is installed
if ! command_exists julia; then
    # Try common locations
    if [ -f "$HOME/.juliaup/bin/julia" ]; then
        export PATH="$HOME/.juliaup/bin:$PATH"
    else
        echo "Julia is not installed. Running setup first..."
        echo ""
        ./setup.sh

        # Re-check after setup
        if [ -f "$HOME/.juliaup/bin/julia" ]; then
            export PATH="$HOME/.juliaup/bin:$PATH"
        fi
    fi
fi

# Launch GUI
echo "Starting SLAM.jl GUI..."
julia --project=. scripts/run_gui.jl
