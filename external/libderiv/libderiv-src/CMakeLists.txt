cmake_minimum_required(VERSION 3.0)
project(libderiv LANGUAGES C CXX)
find_package(libint REQUIRED COMPONENTS ${LIBINT_OPT_AM})
add_subdirectory(libderiv)
