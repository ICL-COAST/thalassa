# Enable FetchContent
include(FetchContent)

# # Add SOFA
# Find an existing SOFA installation
find_package(sofa-f QUIET)
# Download SOFA if an existing installation is not found
if(NOT ${sofa-f_FOUND})
    # Disable the ANSI C version of SOFA
    set(SOFA_BUILD_C OFF)

    # Fetch SOFA
    FetchContent_Declare(sofa
        GIT_REPOSITORY https://github.com/maxhlc/sofa.git
        GIT_TAG 2021-05-12-rev2
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
    FetchContent_Declare(spice
        GIT_REPOSITORY https://github.com/maxhlc/spice.git
        GIT_TAG N0067-rev1
    )

    # Make SPICE available
    FetchContent_MakeAvailable(spice)
endif(NOT ${spice_FOUND})