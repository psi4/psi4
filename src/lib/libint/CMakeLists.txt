#
# Figure out which files libint will produce, and compile them
#
math(EXPR LIBINT_NEW_AM ${LIBINT_OPT_AM}*2)
configure_file(${CMAKE_SOURCE_DIR}/src/lib/libint/libint_config.h.in ${CMAKE_BINARY_DIR}/src/lib/libint/libint_config.h)

include(SetupLibint)
build_libint_files(LI_SRC ${LIBINT_OPT_AM})
set(LI_SRC ${LI_SRC} init_libint.cc vrr_build.c)

#foreach (f ${LI_SRC})
#    message(STATUS "LIBINT CONTAINS ${f}")
#endforeach()
set_libint_compiler_flags()
add_library(int STATIC ${LI_SRC})

#
# Build the compiler
#
set(COMPILER_SRC build_libint.c emit_hrr_build.c emit_hrr_build_macro.c emit_order.c
                 emit_vrr_build.c emit_vrr_build_macro.c mem_man.c constants.c)
add_executable(libint_compiler ${COMPILER_SRC})

#
# Make scratch
#
#add_custom_command(
#    COMMAND ${CMAKE_COMMAND} -E make_directory ${LIBINT_TMPDIR}
#    OUTPUT ${LIBINT_TMPDIR}
#)


# Generate the libint files, after building the compiler, only if libint.compiled doesn't exist
add_custom_command(
    COMMAND ${CMAKE_COMMAND} -E copy ${PROJECT_SOURCE_DIR}/src/lib/libint/vrr_build.c ${PROJECT_BINARY_DIR}/src/lib/libint
    COMMAND libint_compiler
    COMMAND touch libint.compiled
    DEPENDS libint_compiler
    OUTPUT libint.compiled
    OUTPUT ${LI_SRC}
    COMMENT "Generating libint files"
)

# This target is always executed
add_custom_target(libint ALL DEPENDS libint.compiled)
set_property(GLOBAL APPEND PROPERTY LIBLIST int)
install_list_FILES("${PROJECT_BINARY_DIR}/src/lib/libint/libint.h" include/libint)
