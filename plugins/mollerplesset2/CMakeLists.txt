set(CCSD "${CMAKE_CURRENT_SOURCE_DIR}")
set(CCBD "${CMAKE_CURRENT_BINARY_DIR}")
set(PSIEXE ${PROJECT_BINARY_DIR}/bin/psi4${CMAKE_EXECUTABLE_SUFFIX})

set(PLUG "mollerplesset2")
set(PLUGSRC __init__.py globals.h mp2.cc mp2_restricted.cc mp2_unrestricted.cc pymodule.py)
#inputalt.dat
#input.dat
#Makefile
#doc.rst

foreach(src ${PLUGSRC})
    execute_process(COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CCSD}/${src} ${CCBD}/${src})
endforeach(src)

add_custom_target(plugin_${PLUG}
    ALL
    DEPENDS ${PSIEXE}
    COMMAND ${PSIEXE} --new-plugin-makefile
    COMMAND make
    COMMENT "Build ${PLUG} example plugin"
    VERBATIM)

