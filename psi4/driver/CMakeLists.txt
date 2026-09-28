add_subdirectory(p4util)
add_subdirectory(procrouting)
# add_subdirectory(qcdb)
# add_subdirectory(util)

list(APPEND sources
  dependency_check
  driver_findif
  driver_nbody_helper
  gaussian_n
  molutil
  pluginutil
  psifiles
  wrapper_autofrag
  endorsed_plugins
  aliases
  diatomic
  driver_cbs
  driver_nbody
  driver_util
  frac
  inputparser
  qmmm
  wrapper_database
  driver
  schema_wrapper
)

foreach(py ${sources})
    execute_process(COMMAND ${Python_EXECUTABLE} -c
        "from Cython.Build import cythonize; cythonize('${CMAKE_CURRENT_SOURCE_DIR}/${py}.py', language_level=3)"
        RESULT_VARIABLE _exec_python_err)
    # message(${_exec_python_err})
    file(COPY ${py}.c DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
    file(REMOVE ${py}.c)

    add_library(${py} SHARED "${py}.c")
    target_include_directories(${py} PRIVATE "${Python_INCLUDE_DIRS}")
    set_target_properties(${py} PROPERTIES PREFIX "${PYTHON_MODULE_PREFIX}"  # for python module
                                           OUTPUT_NAME ${py}
                                           SUFFIX "${PYTHON_MODULE_EXTENSION}")  # for python module

    install(TARGETS ${py}
            LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}${PYMOD_INSTALL_LIBDIR}/psi4/driver)
endforeach()

# not going to cythonize qcdb/ or util/ . inits would need a special naming procedure or cmake's confused so copying them plaintext
install(DIRECTORY .
        DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}${PYMOD_INSTALL_LIBDIR}/psi4/driver
        FILES_MATCHING PATTERN "__init__.py")

install(DIRECTORY qcdb
        DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}${PYMOD_INSTALL_LIBDIR}/psi4/driver
        FILES_MATCHING PATTERN "*.py")

install(DIRECTORY util
        DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}${PYMOD_INSTALL_LIBDIR}/psi4/driver
        FILES_MATCHING PATTERN "*.py")
