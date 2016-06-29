cmake_minimum_required(VERSION 3.0)
project(dkh LANGUAGES Fortran)

add_library(dkh dkh2-dk4_main.F90)
install(TARGET dkh DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
