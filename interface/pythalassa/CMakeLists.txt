# Enable FetchContent
include(FetchContent)

## Add pybind11
# Find existing pybind11 installation
find_package(pybind11 QUIET)
# Download pybind11 if an existing installation is not found
if(NOT ${pybind11_FOUND})
    # Fetch pybind11
    FetchContent_Declare(pybind11
        GIT_REPOSITORY https://github.com/pybind/pybind11.git
        GIT_TAG v2.10.4
    )

    # Make pybind11 available
    FetchContent_MakeAvailable(pybind11)
endif(NOT ${pybind11_FOUND})

# Collect source files
file(GLOB_RECURSE src_cpp CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp)

# Build Python interface
pybind11_add_module(pythalassa ${src_cpp})
target_link_libraries(pythalassa PUBLIC cthalassa)