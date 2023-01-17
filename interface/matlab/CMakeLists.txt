# Find Matlab
find_package(Matlab)

# Build MEX if Matlab found
if(Matlab_FOUND)
    matlab_add_mex(NAME mthalassa
                   SRC ${CMAKE_CURRENT_SOURCE_DIR}/src/mthalassa.c
                   LINK_TO cthalassa
                   DOCUMENTATION ${CMAKE_CURRENT_SOURCE_DIR}/src/mthalassa.txt
                   R2018a
    )
endif()