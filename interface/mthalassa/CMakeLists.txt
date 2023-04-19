# Find Matlab
find_package(Matlab)

# Warning message for Apple devices
if(APPLE)
    message(WARNING "Issues may occur when compiling MTHALASSA for Matlab on MacOS. Please refer to the documentation.")
endif()

# Collect header files
set(src_headers ${CMAKE_CURRENT_SOURCE_DIR}/include)
file(GLOB_RECURSE src_h CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/include/*.h)

# Collect source files
file(GLOB_RECURSE src_c CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/src/*.c)

if(Matlab_FOUND)
    # Build MEX if Matlab found
    matlab_add_mex(NAME mthalassa
                   SRC ${src_c}
                   LINK_TO cthalassa
                   DOCUMENTATION ${CMAKE_CURRENT_SOURCE_DIR}/mthalassa.txt
                   R2018a
    )

    # Include header files
    target_include_directories(mthalassa PUBLIC ${src_headers})
endif()