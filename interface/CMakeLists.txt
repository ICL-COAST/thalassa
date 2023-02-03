# Add CTHALASSA
if(THALASSA_BUILD_CTHALASSA)
    add_subdirectory(cpp)
endif()

# Add MTHALASSA
if(THALASSA_BUILD_MTHALASSA)
    add_subdirectory(matlab)
endif()