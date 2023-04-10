# Set minimum CMake version
cmake_minimum_required(VERSION 3.16.3)

# Set project metadata
project(thalassa LANGUAGES C CXX Fortran)

# Set C compiler flags
set(CMAKE_C_FLAGS                "")
set(CMAKE_C_FLAGS_DEBUG          "-g -Wall")
set(CMAKE_C_FLAGS_RELEASE        "-O3 -march=native")
set(CMAKE_C_FLAGS_RELWITHDEBINFO "${CMAKE_Fortran_FLAGS_DEBUG} ${CMAKE_CXX_FLAGS_RELEASE}")
set(CMAKE_C_FLAGS_MINSIZEREL     "")

# Set C++ compiler flags
set(CMAKE_CXX_FLAGS                "")
set(CMAKE_CXX_FLAGS_DEBUG          "-g -Wall")
set(CMAKE_CXX_FLAGS_RELEASE        "-O3 -march=native")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_Fortran_FLAGS_DEBUG} ${CMAKE_CXX_FLAGS_RELEASE}")
set(CMAKE_CXX_FLAGS_MINSIZEREL     "")

# Set Fortran compiler flags
set(CMAKE_Fortran_FLAGS                "")
set(CMAKE_Fortran_FLAGS_DEBUG          "-g -Wall")
set(CMAKE_Fortran_FLAGS_RELEASE        "-march=native")
set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "${CMAKE_Fortran_FLAGS_DEBUG} ${CMAKE_Fortran_FLAGS_RELEASE}")
set(CMAKE_Fortran_FLAGS_MINSIZEREL     "")

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Set output directories
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

# Set directory for .mod files
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod)

# Declare options
option(THALASSA_BUILD_CTHALASSA "Build CTHALASSA interface" ON)
option(THALASSA_BUILD_MTHALASSA "Build MTHALASSA interface" ON)

# Force CTHALASSA build if MTHALASSA requested
if(THALASSA_BUILD_MTHALASSA AND (NOT THALASSA_BUILD_CTHALASSA))
    message(WARNING "Forcing the build of CTHALASSA as it is required for MTHALASSA")
    set(THALASSA_BUILD_CTHALASSA ON)
endif()

# Add external libraries
add_subdirectory(external)

# Compile THALASSA library
add_subdirectory(src)

# Compile THALASSA interfaces
add_subdirectory(interface)

# Compile THALASSA applications
add_subdirectory(app)
