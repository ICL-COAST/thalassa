include(FetchContent)

# Show download progress
set(FETCHCONTENT_QUIET FALSE)

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
    spice-fortran
    URL ${SPICE_URL}
)
FetchContent_MakeAvailable(spice-fortran)

# Extract filepaths for source files
FetchContent_GetProperties(spice-fortran)
file(GLOB_RECURSE src_f ${spice-fortran_SOURCE_DIR}/src/spicelib/*.f)

# Set compiler flags
set(CMAKE_Fortran_FLAGS                "-m64 -fPIC -c -std=legacy -Wno-character-truncation")
set(CMAKE_Fortran_FLAGS_DEBUG          "-g")
set(CMAKE_Fortran_FLAGS_RELEASE        "")
set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "${CMAKE_Fortran_FLAGS_DEBUG} ${CMAKE_Fortran_FLAGS_RELEASE}")
set(CMAKE_Fortran_FLAGS_MINSIZEREL     "")

# Build SPICE library
add_library(spice-fortran STATIC ${src_f})
set_target_properties(spice-fortran PROPERTIES OUTPUT_NAME spice)