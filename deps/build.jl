# Build script for SLAM.jl Fortran dependencies
# Compiles coul90 and wrapper into a shared library

using Libdl

# Get the directory containing this script
const depsdir = @__DIR__
const srcdir = depsdir

# Output library name
const libname = Sys.iswindows() ? "libcoul90.dll" :
                Sys.isapple() ? "libcoul90.dylib" : "libcoul90.so"

const libpath = joinpath(depsdir, libname)

# Fortran compiler (try gfortran first)
const FC = get(ENV, "FC", "gfortran")

# Compiler flags
const FFLAGS = ["-O3", "-fPIC", "-shared"]

# Source files (order matters: coul90.f must come before wrapper)
const sources = [
    joinpath(srcdir, "coul90.f"),
    joinpath(srcdir, "coul90_wrapper.f90")
]

function build()
    println("Building Fortran Coulomb library...")
    println("  Compiler: $FC")
    println("  Output: $libpath")

    # Check if source files exist
    for src in sources
        if !isfile(src)
            error("Source file not found: $src")
        end
    end

    # Compile command
    cmd = `$FC $FFLAGS -o $libpath $sources`
    println("  Command: $cmd")

    try
        run(cmd)
        println("Build successful!")
    catch e
        @warn "Build failed: $e"
        @warn "You may need to install gfortran or set the FC environment variable"
        rethrow()
    end

    # Verify the library was created
    if isfile(libpath)
        println("Library created: $libpath")
        # Try to load it
        try
            handle = Libdl.dlopen(libpath)
            Libdl.dlclose(handle)
            println("Library loads successfully!")
        catch e
            @warn "Library created but failed to load: $e"
        end
    else
        error("Library was not created")
    end
end

# Run build if this script is executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    build()
end
