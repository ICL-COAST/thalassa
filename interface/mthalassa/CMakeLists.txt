# Find Matlab
find_package(Matlab)

# Warning message for Apple devices
if(APPLE)
    message(WARNING "Issues may occur when compiling MTHALASSA for Matlab on MacOS. Please refer to the documentation.")
endif()

# Build MEX if Matlab found
if(Matlab_FOUND)
    matlab_add_mex(NAME mthalassa
                   SRC ${CMAKE_CURRENT_SOURCE_DIR}/src/mthalassa.c
                   LINK_TO cthalassa
                   DOCUMENTATION ${CMAKE_CURRENT_SOURCE_DIR}/mthalassa.txt
                   R2018a
    )
endif()