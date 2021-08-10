add_subdirectory(dft)
add_subdirectory(findif_response_utils)
# add_subdirectory(libcubeprop)
# add_subdirectory(mcscf)
add_subdirectory(response)
add_subdirectory(sapt)
add_subdirectory(scf_proc)
add_subdirectory(solvent)

list(APPEND sources
  empirical_dispersion
  interface_cfour
  proc_table
  roa
  wrappers_cfour
  proc_util
  proc
)

foreach(py ${sources})
    execute_process(COMMAND ${Python_EXECUTABLE} -c
        "from Cython.Build import cythonize; cythonize('${CMAKE_CURRENT_SOURCE_DIR}/${py}.py', language_level=3)"
        RESULT_VARIABLE _exec_python_err)
    file(COPY ${py}.c DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
    file(REMOVE ${py}.c)

    add_library(${py} SHARED "${py}.c")
    target_include_directories(${py} PRIVATE "${Python_INCLUDE_DIRS}")
    set_target_properties(${py} PROPERTIES PREFIX "${PYTHON_MODULE_PREFIX}"  # for python module
                                           OUTPUT_NAME ${py}
                                           SUFFIX "${PYTHON_MODULE_EXTENSION}")  # for python module
 
    install(TARGETS ${py}
            LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}${PYMOD_INSTALL_LIBDIR}/psi4/driver/procrouting)
endforeach()

install(DIRECTORY libcubeprop
        DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}${PYMOD_INSTALL_LIBDIR}/psi4/driver/procrouting
        FILES_MATCHING PATTERN "*.py")

install(DIRECTORY mcscf
        DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}${PYMOD_INSTALL_LIBDIR}/psi4/driver/procrouting
        FILES_MATCHING PATTERN "*.py")
