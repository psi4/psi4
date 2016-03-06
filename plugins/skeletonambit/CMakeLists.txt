# check out and build fresh plugin
# creates a dir skeleton w/i skeleton then copies files up to satisfy requirements of both CMake and Psi4 plugins
# symlinks input.dat back to source dir to satisfy TestingMacro
set(CCSD "${CMAKE_CURRENT_SOURCE_DIR}")
set(CCBD "${CMAKE_CURRENT_BINARY_DIR}")
set(PSIEXE ${PROJECT_BINARY_DIR}/bin/psi4${CMAKE_EXECUTABLE_SUFFIX})

set(PLUG "skeletonambit")
set(PLUGSRC ${PLUG}.cc input.dat Makefile pymodule.py __init__.py doc.rst)

add_custom_target(plugin_${PLUG}
    ALL
    DEPENDS ${PSIEXE}
    COMMAND ${CMAKE_COMMAND} -E remove_directory ${CCBD}/${PLUG}
    COMMAND ${PSIEXE} -l ${PROJECT_SOURCE_DIR}/share --new-plugin ${PLUG}
    COMMAND ${CMAKE_COMMAND} -E copy_directory ${CCBD}/${PLUG} ${CCBD}
    COMMAND ${CMAKE_COMMAND} -E remove_directory ${CCBD}/${PLUG}
    COMMAND $(MAKE)
    COMMAND ${CMAKE_COMMAND} -E create_symlink ${CCBD}/input.dat ${CCSD}/input.dat
    COMMENT "Build ${PLUG} example plugin"
    VERBATIM)

include(TestingMacros)
add_regression_test(${PLUG} "psi;plug")

