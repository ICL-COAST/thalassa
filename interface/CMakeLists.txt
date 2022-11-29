# Enable C
enable_language(C)

# Add CTHALASSA executable and link to THALASSA
add_executable(cthalassa cthalassa.c cthalassa.f90)
target_link_libraries(cthalassa PUBLIC thalassa)