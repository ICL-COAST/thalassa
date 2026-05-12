# Add Eigen
# TODO: scope build testing
set(BUILD_TESTING OFF)
set(EIGEN_BUILD_TESTING OFF)
set(EIGEN_MPL2_ONLY ON)
set(EIGEN_BUILD_PKGCONFIG OFF)
set(EIGEN_BUILD_DOC OFF)
add_subdirectory(eigen)

# Add SOFA
set(SOFA_BUILD_C OFF)
add_subdirectory(sofa)

# Add SPICE
add_subdirectory(spice)
