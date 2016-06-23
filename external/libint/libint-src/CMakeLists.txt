cmake_minimum_required(VERSION 3.0)
project(libint LANGUAGES C CXX)
add_subdirectory(libint)
configure_file("libintConfig.cmake.in" "libintConfig.cmake" @ONLY)
install(FILES ${CMAKE_BINARY_DIR}/libintConfig.cmake 
        DESTINATION ${CMAKE_INSTALL_PREFIX})
