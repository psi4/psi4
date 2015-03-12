include(FindGit)

if(GIT_FOUND)
   add_custom_target(
        git_update
        COMMAND ${GIT_EXECUTABLE} submodule init
        COMMAND ${GIT_EXECUTABLE} submodule update
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})
else()
    if(DEVELOPMENT_CODE)
       message("-- Git not found. You need Git for the Git submodule mechanism to work.")
    endif()
endif()

include(ExternalProject)

macro(add_external _project _external_project_cmake_args _testing_command)
    if(DEVELOPMENT_CODE AND GIT_FOUND)
       set(UPDATE_COMMAND ${GIT_EXECUTABLE} submodule update)
    else()
       set(UPDATE_COMMAND echo)
    endif()

    add_custom_target(
        check_external_timestamp_${_project}
	COMMAND ${PYTHON_EXECUTABLE} ${PROJECT_SOURCE_DIR}/cmake/check_external_timestamp.py
                       ${PROJECT_BINARY_DIR}/interfaces/${_project}-stamp/${_project}-configure
                       ${PROJECT_BINARY_DIR}/interfaces/${_project}-stamp
                       ${PROJECT_SOURCE_DIR}/interfaces/${_project}
    )

    if(_testing_command)
       ExternalProject_Add(${_project}
           DOWNLOAD_COMMAND ${UPDATE_COMMAND}
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
    else()
       ExternalProject_Add(${_project}
           DOWNLOAD_COMMAND ${UPDATE_COMMAND}
           DOWNLOAD_DIR ${PROJECT_SOURCE_DIR}
           SOURCE_DIR ${PROJECT_SOURCE_DIR}/interfaces/${_project}
           BINARY_DIR ${PROJECT_BINARY_DIR}/interfaces/${_project}-build
           STAMP_DIR ${PROJECT_BINARY_DIR}/interfaces/${_project}-stamp
           TMP_DIR ${PROJECT_BINARY_DIR}/interfaces/${_project}-tmp
           INSTALL_DIR ${PROJECT_BINARY_DIR}/interfaces
           CMAKE_ARGS ${_external_project_cmake_args}
           )
    endif()

    link_directories(${PROJECT_BINARY_DIR}/interfaces/lib)
    link_directories(${PROJECT_BINARY_DIR}/interfaces/${_project}-build/interfaces/lib)

    if(DEVELOPMENT_CODE)
       add_dependencies(${_project} git_update)
       add_dependencies(check_external_timestamp_${_project} git_update)
    endif()
    add_dependencies(${_project} check_external_timestamp_${_project})
endmacro()
