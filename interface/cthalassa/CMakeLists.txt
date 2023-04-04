# Define CTHALASSA options
option(CTHALASSA_USE_FORK "Use fork in CTHALASSA to spawn subprocesses" ON)

# Disable fork on non-UNIX systems
if(NOT UNIX)
    set(CTHALASSA_USE_FORK OFF)
endif()

# Collect header files
set(src_headers ${CMAKE_CURRENT_SOURCE_DIR}/include)
file(GLOB_RECURSE src_h CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/include/*.h)
file(GLOB_RECURSE src_hpp CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/include/*.hpp)

# Collect source files
file(GLOB_RECURSE src_c CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/src/*.c)
file(GLOB_RECURSE src_cpp CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp)
file(GLOB_RECURSE src_f90 CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/src/*.f90)

# Add CTHALASSA library
add_library(cthalassa STATIC ${src_c} ${src_cpp} ${src_h} ${src_hpp} ${src_f90})

# Include CTHALASSA headers
target_include_directories(cthalassa PUBLIC ${src_headers})

# Link CTHALASSA to THALASSA
target_link_libraries(cthalassa PUBLIC thalassa)

# Enable threading
find_package(Threads REQUIRED)
target_link_libraries(cthalassa PUBLIC Threads::Threads)

# Enable fork for propagations
if(CTHALASSA_USE_FORK)
    target_compile_definitions(cthalassa PUBLIC "CTHALASSA_USE_FORK")
endif()

# Add CTHALASSA interface example
add_executable(cthalassa_example cthalassa_example.cpp)

# Link to CTHALASSA
target_link_libraries(cthalassa_example PRIVATE cthalassa)