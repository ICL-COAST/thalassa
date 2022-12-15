# Enable C
enable_language(C)

# Collect files
file(GLOB_RECURSE src_c CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/src/*.c)
file(GLOB_RECURSE src_h CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/include/*.h)
file(GLOB_RECURSE src_hpp CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/include/*.hpp)
file(GLOB_RECURSE src_f90 CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/src/*.f90)
set(src_headers ${CMAKE_CURRENT_SOURCE_DIR}/include)

# Add CTHALASSA executable, include headers, and link to THALASSA
add_library(cthalassa ${src_c} ${src_h} ${src_hpp} ${src_f90})
target_include_directories(cthalassa PUBLIC ${src_headers})
target_link_libraries(cthalassa PUBLIC thalassa)