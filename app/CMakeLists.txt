# Compile THALASSA main file
add_executable(${PROJECT_NAME}_main ${CMAKE_CURRENT_SOURCE_DIR}/thalassa_main.f90)
target_link_libraries(${PROJECT_NAME}_main ${PROJECT_NAME})