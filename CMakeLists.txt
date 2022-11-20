# Set minimum CMake version
cmake_minimum_required(VERSION 3.16.3)

# Set project metadata
project(thalassa)

# Enable Fortran
enable_language(Fortran)

# Set output directories
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

# TODO: add option to download SPICE (w/ kernels) and SOFA

# Add SPICE and SOFA
add_library(sofa STATIC IMPORTED)
set_property(TARGET sofa PROPERTY IMPORTED_LOCATION ${CMAKE_CURRENT_SOURCE_DIR}/external/libsofa.a)
add_library(spice STATIC IMPORTED)
set_property(TARGET spice PROPERTY IMPORTED_LOCATION ${CMAKE_CURRENT_SOURCE_DIR}/external/spicelib.a) 

# Collect filepaths
file(GLOB_RECURSE src_f ${CMAKE_CURRENT_SOURCE_DIR}/src/*.f)
file(GLOB_RECURSE src_for ${CMAKE_CURRENT_SOURCE_DIR}/src/*.for)
file(GLOB_RECURSE src_f90 ${CMAKE_CURRENT_SOURCE_DIR}/src/*.f90)

# Set file-specific compiler flags
set_source_files_properties(${src_f} PROPERTIES COMPILE_FLAGS "-fdefault-real-8 -std=legacy")
set_source_files_properties(${src_for} PROPERTIES COMPILE_FLAGS "-fdefault-real-8 -std=legacy -w")

# Compile THALASSA library
add_library(${PROJECT_NAME} STATIC ${src_f} ${src_for} ${src_f90})
target_link_libraries(${PROJECT_NAME} PUBLIC sofa spice)

# Compile THALASSA main file
add_executable(${PROJECT_NAME}_main ${CMAKE_CURRENT_SOURCE_DIR}/src/thalassa_main.f90)
target_link_libraries(${PROJECT_NAME}_main ${PROJECT_NAME})
