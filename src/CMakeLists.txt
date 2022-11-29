# Collect filepaths
file(GLOB_RECURSE src_f ${CMAKE_CURRENT_SOURCE_DIR}/*.f)
file(GLOB_RECURSE src_for ${CMAKE_CURRENT_SOURCE_DIR}/*.for)
file(GLOB_RECURSE src_f90 ${CMAKE_CURRENT_SOURCE_DIR}/*.f90)

# Set file-specific compiler flags
set_source_files_properties(${src_f} PROPERTIES COMPILE_FLAGS "-fdefault-real-8 -std=legacy")
set_source_files_properties(${src_for} PROPERTIES COMPILE_FLAGS "-fdefault-real-8 -std=legacy -w")

# Compile THALASSA library
add_library(${PROJECT_NAME} STATIC ${src_f} ${src_for} ${src_f90})
target_link_libraries(${PROJECT_NAME} PUBLIC sofa spice)