include(FetchContent)

# Show download progress
set(FETCHCONTENT_QUIET FALSE)

# Set URL based on operating system
if(UNIX)
    # Version with UNIX line terminations
    set(SOFA_URL http://www.iausofa.org/2021_0512_F/sofa_f-20210512.tar.gz)
elseif(WIN32)
    # Version with DOS line terminations
    set(SOFA_URL https://www.iausofa.org/2021_0512_F/sofa_f-20210512.zip)
endif()

# Download SOFA source
FetchContent_Declare(
    sofa-fortran
    URL ${SOFA_URL}
)
FetchContent_MakeAvailable(sofa-fortran)

# Extract filepaths for source files
FetchContent_GetProperties(sofa-fortran)
file(GLOB_RECURSE src_for ${sofa-fortran_SOURCE_DIR}/20210512/f77/src/*.for)

# Set file-specific compiler flags
set_source_files_properties(${src_for} PROPERTIES COMPILE_FLAGS "")

# Build SOFA library
add_library(sofa-fortran STATIC ${src_for})
set_target_properties(sofa-fortran PROPERTIES OUTPUT_NAME sofa)