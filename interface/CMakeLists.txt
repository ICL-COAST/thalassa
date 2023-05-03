# Add CTHALASSA
if(THALASSA_BUILD_CTHALASSA)
    add_subdirectory(cthalassa)
endif()

# Add MTHALASSA
if(THALASSA_BUILD_MTHALASSA)
    add_subdirectory(mthalassa)
endif()

# Add PYTHALASSA
if(THALASSA_BUILD_PYTHALASSA)
    add_subdirectory(pythalassa)
endif()