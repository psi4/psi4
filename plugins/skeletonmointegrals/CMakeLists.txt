# check out and build fresh plugin
# creates a dir skeleton w/i skeleton then copies files up to satisfy requirements of both CMake and Psi4 plugins
# symlinks input.dat back to source dir to satisfy TestingMacro
set(CCSD "${CMAKE_CURRENT_SOURCE_DIR}")
set(CCBD "${CMAKE_CURRENT_BINARY_DIR}")
set(PSIEXE ${STAGED_INSTALL_PREFIX}/${CMAKE_INSTALL_BINDIR}/psi4)

set(PLUG "skeletonmointegrals")

add_custom_target(plugin_${PLUG}
    ALL
    DEPENDS psi4-core
    COMMAND ${CMAKE_COMMAND} -E remove_directory ${CCBD}/${PLUG}
    COMMAND ${PSIEXE} --plugin-name ${PLUG} --plugin-template mointegrals
    COMMAND ${CMAKE_COMMAND} -E chdir "${CCBD}/${PLUG}" cmake -C ${STAGED_INSTALL_PREFIX}/share/cmake/psi4/psi4PluginCache.cmake -DCMAKE_PREFIX_PATH=${STAGED_INSTALL_PREFIX} .
    COMMAND ${CMAKE_COMMAND} -E chdir "${CCBD}/${PLUG}" ${CMAKE_MAKE_PROGRAM}
    COMMAND ${CMAKE_COMMAND} -E create_symlink ${CCBD}/${PLUG}/input.dat ${CCSD}/input.dat
    COMMAND ${CMAKE_COMMAND} -E create_symlink "${PLUG}/${PLUG}.so" "${PLUG}.so"
    COMMAND ${CMAKE_COMMAND} -E create_symlink "${PLUG}/__init__.py" "__init__.py"
    COMMAND ${CMAKE_COMMAND} -E create_symlink "${PLUG}/pymodule.py" "pymodule.py"
    COMMENT "Build ${PLUG} example plugin"
    VERBATIM)

include(TestingMacros)
add_regression_test(${PLUG} "psi;plug")

