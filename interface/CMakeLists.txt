# Enable C
enable_language(C)

# Collect files
file(GLOB_RECURSE src_c ${CMAKE_CURRENT_SOURCE_DIR}/src/*.c)
file(GLOB_RECURSE src_f90 ${CMAKE_CURRENT_SOURCE_DIR}/src/*.f90)
set(src_h ${CMAKE_CURRENT_SOURCE_DIR}/include)

# Add CTHALASSA executable, include headers, and link to THALASSA
add_executable(cthalassa ${src_c} ${src_f90})
target_include_directories(cthalassa PUBLIC ${src_h})
target_link_libraries(cthalassa PUBLIC thalassa)