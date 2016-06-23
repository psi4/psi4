cmake_minimum_required(VERSION 3.0)
project(libderiv_compiler LANGUAGES C)
find_package(libint REQUIRED COMPONENTS ${LIBINT_OPT_AM})
include_directories(${LIBINT_INCLUDE_DIRS})
list(APPEND CMAKE_C_FLAGS -lm)
math(EXPR LIBDERIV_OPT_AM1 ${LIBINT_OPT_AM}-1) # A.M. level for 1st derivative ERIs
math(EXPR LIBDERIV_OPT_AM2 ${LIBINT_OPT_AM}-2) # A.M. level for 2nd derivative ERIs
math(EXPR LIBDERIV_NEW_AM1 ${LIBDERIV_OPT_AM1}*2)
math(EXPR LIBDERIV_NEW_AM2 ${LIBDERIV_OPT_AM2}*2)
math(EXPR LIBDERIV_NEW_AM12 ${LIBDERIV_OPT_AM2}*2)
add_definitions(-DLIBDERIV_NEW_AM1=${LIBDERIV_NEW_AM1})
add_definitions(-DLIBDERIV_NEW_AM2=${LIBDERIV_NEW_AM2})
add_definitions(-DLIBDERIV_NEW_AM12=${LIBDERIV_NEW_AM2})
#
# Build the compiler
#
set(COMPILER_SRC build_libderiv.c 
                 emit_d1hrr_build_macro.c 
                 emit_deriv1_managers.c 
                 emit_deriv_build_macro.c
                 emit_d1hrr_build.c
                 emit_deriv12_managers.c
                 emit_deriv_build.c
                 mem_man.c
)
add_executable(libderiv_compiler ${COMPILER_SRC})
install(TARGETS libderiv_compiler DESTINATION ${CMAKE_INSTALL_PREFIX})
