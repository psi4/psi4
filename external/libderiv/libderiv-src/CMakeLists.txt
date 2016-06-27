cmake_minimum_required(VERSION 3.0)
project(libderiv LANGUAGES C CXX)
find_package(libint REQUIRED COMPONENTS ${LIBINT_OPT_AM})
add_subdirectory(libderiv)
configure_file("libderivConfig.cmake.in" "libderivConfig.cmake" @ONLY)
install(FILES ${CMAKE_BINARY_DIR}/libderivConfig.cmake 
        DESTINATION ${CMAKE_INSTALL_PREFIX})