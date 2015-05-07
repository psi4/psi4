set(CCSD "${CMAKE_CURRENT_SOURCE_DIR}")
set(CCBD "${CMAKE_CURRENT_BINARY_DIR}")
set(PSIEXE ${PROJECT_BINARY_DIR}/bin/psi4${CMAKE_EXECUTABLE_SUFFIX})

set(PLUG "backtrans")
set(PLUGSRC backtrans.cc input.dat)

foreach(src ${PLUGSRC})
    execute_process(COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CCSD}/${src} ${CCBD}/${src})
endforeach(src)

add_custom_target(plugin_${PLUG}
    ALL
    DEPENDS ${PSIEXE}
    COMMAND ${PSIEXE} -l ${PROJECT_SOURCE_DIR}/lib --new-plugin-makefile
    COMMAND $(MAKE)
    COMMENT "Build ${PLUG} example plugin"
    VERBATIM)

include(TestingMacros)
add_regression_test(${PLUG} "psi;plug")

