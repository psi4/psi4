add_subdirectory(aointegrals)
add_subdirectory(backtrans)
add_subdirectory(mointegrals)
add_subdirectory(mollerplesset2)
add_subdirectory(sointegrals)

# check out and build fresh plugin
set(PSIEXE ${PROJECT_BINARY_DIR}/bin/psi4${CMAKE_EXECUTABLE_SUFFIX})
set(PLUG "skeleton")
add_custom_target(plugin_${PLUG}
    ALL
    DEPENDS ${PSIEXE}
    COMMAND ${CMAKE_COMMAND} -E remove_directory ${CMAKE_CURRENT_BINARY_DIR}/${PLUG}
    COMMAND ${PSIEXE} --psidatadir ${CMAKE_SOURCE_DIR}/lib --new-plugin ${PLUG}
    COMMAND ${CMAKE_COMMAND} -E chdir ${CMAKE_CURRENT_BINARY_DIR}/${PLUG} $(MAKE)
    COMMENT "Build ${PLUG} example plugin"
    VERBATIM)

include(TestingMacros)
add_regression_test(${PLUG} "psi;plug")

