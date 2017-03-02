# CMake driving pytest copied from pybind11:
#   https://github.com/pybind/pybind11/blob/master/tests/CMakeLists.txt

# Make sure pytest is found or produce a fatal error
if(NOT psi4_PYTEST_FOUND)
    execute_process(COMMAND ${PYTHON_EXECUTABLE} -c "import pytest"
                    OUTPUT_QUIET ERROR_QUIET
                    RESULT_VARIABLE _exec_python_err)
    if(_exec_python_err)
        message(STATUS "Running some tests requires pytest. Please install it manually (try: ${PYTHON_EXECUTABLE} -m pip install pytest)")
    endif()
    set(psi4_PYTEST_FOUND TRUE CACHE INTERNAL "")
endif()

# A single command to remake and run the pytest smoke tests
add_custom_target(pytest
                  COMMAND PYTHONPATH=${STAGED_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}${PYMOD_INSTALL_LIBDIR}
                          ${PYTHON_EXECUTABLE} -m pytest -rws -v --capture=sys ${CMAKE_CURRENT_SOURCE_DIR}
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
