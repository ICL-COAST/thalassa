# Add option for dependencies
option(THALASSA_CONFIGURE_EXTERNAL "Configure external dependencies" OFF)

# Configure SOFA and SPICE
if(THALASSA_CONFIGURE_EXTERNAL)
    execute_process(
        COMMAND sh sofa.sh
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    )
    execute_process(
        COMMAND sh spice.sh
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    )
endif()

# Add SOFA
add_library(sofa STATIC IMPORTED GLOBAL)
set_property(TARGET sofa PROPERTY IMPORTED_LOCATION ${CMAKE_CURRENT_SOURCE_DIR}/libsofa.a)

# Add SPICE
add_library(spice STATIC IMPORTED GLOBAL)
set_property(TARGET spice PROPERTY IMPORTED_LOCATION ${CMAKE_CURRENT_SOURCE_DIR}/spicelib.a)