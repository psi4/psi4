#
# Figure out which files libint will produce, and compile them
#
include(../CMakeLists.libint.txt)
build_libint_files(LI_SRC ${LIBINT_OPT_AM})
set(LI_SRC ${LI_SRC} init_libint.cc vrr_build.c)

#foreach (f ${LI_SRC})
#    message(STATUS "LIBINT CONTAINS ${f}")
#endforeach()
add_library(int STATIC ${LI_SRC})
#
# Build the compiler
#
set(COMPILER_SRC build_libint.c emit_hrr_build.c emit_hrr_build_macro.c emit_order.c
                 emit_vrr_build.c emit_vrr_build_macro.c mem_man.c)
add_executable(libint_compiler ${COMPILER_SRC})

#
# Make scratch
#
#add_custom_command(
#    COMMAND ${CMAKE_COMMAND} -E make_directory ${LIBINT_TMPDIR}
#    OUTPUT ${LIBINT_TMPDIR}
#)


#
# Generate the libint files, after building the compiler, only if libint.compiled doesn't exist
#
add_custom_command(
    COMMAND libint_compiler
    COMMAND touch libint.compiled
    DEPENDS libint_compiler
    OUTPUT libint.compiled
    OUTPUT ${LI_SRC}
)


#
# This target is always executed
#
add_custom_target(
    libint ALL
    DEPENDS libint.compiled
)
