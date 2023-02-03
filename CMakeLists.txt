# Set minimum CMake version
cmake_minimum_required(VERSION 3.16.3)

# Set project metadata
project(thalassa LANGUAGES C CXX Fortran)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)

# Set output directories
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

# Set directory for .mod files
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod)

# Add external libraries
add_subdirectory(external)

# Compile THALASSA library
add_subdirectory(src)

# Compile CTHALASSA library
add_subdirectory(interface)

# Compile THALASSA main file
add_subdirectory(app)
