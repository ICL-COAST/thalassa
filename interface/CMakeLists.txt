# Enable C
enable_language(C)

# Collect files
file(GLOB_RECURSE src_c CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/src/*.c)
file(GLOB_RECURSE src_cpp CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp)
file(GLOB_RECURSE src_h CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/include/*.h)
file(GLOB_RECURSE src_hpp CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/include/*.hpp)
file(GLOB_RECURSE src_f90 CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/src/*.f90)
set(src_headers ${CMAKE_CURRENT_SOURCE_DIR}/include)

# Add CTHALASSA library
add_library(cthalassa ${src_c} ${src_cpp} ${src_h} ${src_hpp} ${src_f90})
# Include CTHALASSA headers
target_include_directories(cthalassa PUBLIC ${src_headers})
# Link CTHALASSA to THALASSA
target_link_libraries(cthalassa PUBLIC thalassa)
# Enable fork for propagations
target_compile_definitions(cthalassa PUBLIC CTHALASSA_USE_FORK)

# Add CTHALASSA interface example
add_executable(interface_example interface_example.cpp)
# Link to CTHALASSA
target_link_libraries(interface_example PRIVATE cthalassa)