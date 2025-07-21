# Enable FetchContent
include(FetchContent)

# Collect filepaths
file(GLOB_RECURSE src_f ${CMAKE_CURRENT_SOURCE_DIR}/*.f)
file(GLOB_RECURSE src_for ${CMAKE_CURRENT_SOURCE_DIR}/*.for)
file(GLOB_RECURSE src_f90 ${CMAKE_CURRENT_SOURCE_DIR}/*.f90)

# Set file-specific compiler flags
set_source_files_properties(${src_f} PROPERTIES COMPILE_FLAGS "-fdefault-real-8 -std=legacy")
set_source_files_properties(${src_for} PROPERTIES COMPILE_FLAGS "-fdefault-real-8 -std=legacy -w")

# Compile THALASSA library
add_library(${PROJECT_NAME} STATIC ${src_f} ${src_for} ${src_f90})

# Ensure PIC
set_target_properties(${PROJECT_NAME} PROPERTIES POSITION_INDEPENDENT_CODE ON)

## Add SOFA
# Find an existing SOFA installation
find_package(sofa-f QUIET)
# Download SOFA if an existing installation is not found
if(NOT ${sofa-f_FOUND})
    # Disable the ANSI C version of SOFA
    set(SOFA_BUILD_C OFF)

    # Fetch SOFA
    FetchContent_Declare(
        sofa
        GIT_REPOSITORY https://github.com/maxhlc/sofa.git
        GIT_TAG 2021-05-12-rev2
        GIT_SHALLOW TRUE
    )

    # Make SOFA available
    FetchContent_MakeAvailable(sofa)
endif(NOT ${sofa-f_FOUND})

## Add SPICE
# Find an existing SPICE installation
find_package(spice QUIET)
# Download SPICE if an existing installation is not found
if(NOT ${spice_FOUND})
    # Fetch SPICE
    FetchContent_Declare(
        spice
        GIT_REPOSITORY https://github.com/maxhlc/spice.git
        GIT_TAG N0067-rev1
        GIT_SHALLOW TRUE
    )

    # Make SPICE available
    FetchContent_MakeAvailable(spice)
endif(NOT ${spice_FOUND})

# Link to dependencies
target_link_libraries(${PROJECT_NAME} PUBLIC sofa-f spice)