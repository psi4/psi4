cmake_minimum_required(VERSION 3.0)
project(libint_compiler LANGUAGES C)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -lm")
math(EXPR LIBINT_NEW_AM ${LIBINT_OPT_AM}*2)
add_definitions(-DLIBINT_NEW_AM=${LIBINT_NEW_AM})
add_definitions(-DLIBINT_OPT_AM=${LIBINT_OPT_AM})
add_definitions(-DLIBINT_MAX_CLASS_SIZE=300)
add_definitions(-DLIBINT_LONG_DOUBLE=0)
#
# Build the compiler
#
set(COMPILER_SRC build_libint.c 
                 emit_hrr_build.c 
                 emit_hrr_build_macro.c 
                 emit_order.c
                 emit_vrr_build.c 
                 emit_vrr_build_macro.c 
                 mem_man.c
                 constants.c
)
add_executable(libint_compiler ${COMPILER_SRC})
install(TARGETS libint_compiler DESTINATION ${CMAKE_INSTALL_PREFIX})
