#
# Build the compiler
#
set(COMPILER_SRC build_libint.c emit_hrr_build.c emit_hrr_build_macro.c emit_order.c
                 emit_vrr_build.c emit_vrr_build_macro.c mem_man.c)# vrr_build.c)
add_executable(libint_compiler ${COMPILER_SRC})

#
# Make scratch
#
add_custom_command(
    COMMAND ${CMAKE_COMMAND} -E make_directory ${LIBINT_TMPDIR}
    OUTPUT ${LIBINT_TMPDIR}
)

#
# Figure out which files libint will produce, and compile them
#
include(../CMakeLists.libint.txt)
build_libint_files(LI_SRC ${LIBINT_OPT_AM} ${LIBINT_TMPDIR})
set(LI_SRC ${INTFILES} ${LIBINT_TMPDIR}/init_libint.cc)
foreach (f ${LI_SRC})
    message(STATUS "${f}")
endforeach()
add_library(libint, ${LI_SRC})

#
# Generate the libint files, after building the compiler and scratch, only if libint.compiled doesn't exist
#
add_custom_command(
    COMMAND libint_compiler
    WORKING_DIRECTORY ${LIBINT_TMPDIR}
    DEPENDS libint_compiler
    DEPENDS ${LIBINT_TMPDIR}
    OUTPUT ${LIBINT_TMPDIR}/init_libint.cc
    OUTPUT libint.compiled
)


#
# This target is always executed
#
add_custom_target(
    libint ALL
    DEPENDS libint.compiled
)
