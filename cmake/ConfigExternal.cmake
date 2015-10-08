include(ExternalProject)

macro(add_external _project _external_project_cmake_args)

    add_custom_target(
        check_external_timestamp_${_project}
        COMMAND ${PYTHON_EXECUTABLE} ${PROJECT_SOURCE_DIR}/cmake/check_external_timestamp.py
                       ${PROJECT_BINARY_DIR}/interfaces/${_project}-stamp/${_project}-configure
                       ${PROJECT_BINARY_DIR}/interfaces/${_project}-stamp
                       ${PROJECT_SOURCE_DIR}/interfaces/${_project}
    )

    if("${_testing_command}" STREQUAL "${ARGV0}")
       ExternalProject_Add(${_project}
           DOWNLOAD_COMMAND echo
           DOWNLOAD_DIR ${PROJECT_SOURCE_DIR}
           SOURCE_DIR ${PROJECT_SOURCE_DIR}/interfaces/${_project}
           BINARY_DIR ${PROJECT_BINARY_DIR}/interfaces/${_project}-build
           STAMP_DIR ${PROJECT_BINARY_DIR}/interfaces/${_project}-stamp
           TMP_DIR ${PROJECT_BINARY_DIR}/interfaces/${_project}-tmp
           INSTALL_DIR ${PROJECT_BINARY_DIR}/interfaces
           CMAKE_ARGS ${_external_project_cmake_args}
           )
    else()
       # For unfathomable reasons, CMake expects the TEST_COMMAND to be ;-separated list...
       separate_arguments(_testing_command)
       ExternalProject_Add(${_project}
           DOWNLOAD_COMMAND echo
           DOWNLOAD_DIR ${PROJECT_SOURCE_DIR}
           SOURCE_DIR ${PROJECT_SOURCE_DIR}/interfaces/${_project}
           BINARY_DIR ${PROJECT_BINARY_DIR}/interfaces/${_project}-build
           STAMP_DIR ${PROJECT_BINARY_DIR}/interfaces/${_project}-stamp
           TMP_DIR ${PROJECT_BINARY_DIR}/interfaces/${_project}-tmp
           INSTALL_DIR ${PROJECT_BINARY_DIR}/interfaces
           CMAKE_ARGS ${_external_project_cmake_args}
           TEST_BEFORE_INSTALL 1
           TEST_COMMAND "${_testing_command}"
           LOG_TEST 1
           )
   endif()

    link_directories(${PROJECT_BINARY_DIR}/interfaces/lib)
    link_directories(${PROJECT_BINARY_DIR}/interfaces/${_project}-build/interfaces/lib)

    add_dependencies(${_project} check_external_timestamp_${_project})

endmacro()
