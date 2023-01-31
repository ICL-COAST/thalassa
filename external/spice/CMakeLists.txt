include(FetchContent)

# Set URL based on operating system
if(UNIX)
    # Set Linux version
    set(SPICE_URL https://naif.jpl.nasa.gov/pub/naif/toolkit//FORTRAN/PC_Linux_gfortran_64bit/packages/toolkit.tar.Z)
elseif(APPLE)
    # Set Apple version
    set(SPICE_URL https://naif.jpl.nasa.gov/pub/naif/toolkit//FORTRAN/MacIntel_OSX_gfortran_64bit/packages/toolkit.tar.Z)
endif()

# Download SPICE source
FetchContent_Declare(
    spice
    URL ${SPICE_URL}
)
FetchContent_MakeAvailable(spice)

# Extract filepaths for source files
FetchContent_GetProperties(spice)
file(GLOB_RECURSE src_f ${spice_SOURCE_DIR}/src/spicelib/*.f)

# Set file-specific compiler flags
set_source_files_properties(${src_f} PROPERTIES COMPILE_FLAGS "-m64 -fPIC -c -std=legacy -Wno-character-truncation")

# Build SPICE library
add_library(spice STATIC ${src_f})