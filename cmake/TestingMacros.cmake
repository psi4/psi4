# This macro is used to add a regression test.
# Tests are identified by their name (one) and by labels (more than one)
# all tests with label quicktests should run in less than 30 seconds sequential on a 2GHz CPU
macro(add_regression_test _name _labels)
    # _labels is not a list, it's a string... Transform it into a list
    set(labels)
    string(REPLACE ";" " " _labels "${_labels}")
    foreach(_label "${_labels}")
        list(APPEND labels ${_label})
    endforeach()
    unset(_labels)
    # This is the psi command to actually run the tests
    set(PSIEXE ${STAGED_INSTALL_PREFIX}/${CMAKE_INSTALL_BINDIR}/psi4)
    # This is the path to the psi4 library
    set(PSILIB ${STAGED_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}${PYMOD_INSTALL_LIBDIR})
    # This is the python script that we call, to call psi4, to run the tests
    set(TESTEXE ${PROJECT_SOURCE_DIR}/tests/runtest.py)
    # This is the psidatadir directory that the script that we call, to call psi4, to run the tests gets its info
    set(PSIDATADIR ${STAGED_INSTALL_PREFIX}/${CMAKE_INSTALL_DATADIR}/psi4)

    # A full report
    set(LOGFILE ${PROJECT_BINARY_DIR}/testresults.log)

    # Generic setup
    set(TEST_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR})
    # Some tests are in subdirectories of subdirectories (cfour, mrcc, dftd3) 
    get_filename_component(dir ${TEST_SRC_DIR} PATH)
    get_filename_component(dir ${dir} NAME)
    if("${dir}" STREQUAL "tests")
        set(TEST_RUN_DIR ${PROJECT_BINARY_DIR}/tests/${_name})
    elseif(("${dir}" STREQUAL "plugins") AND ("${_name}" MATCHES "skeleton"))
        set(TEST_RUN_DIR ${PROJECT_BINARY_DIR}/tests/plugins/${_name})
        set(TEST_SRC_DIR ${TEST_RUN_DIR}/${_name})
    else()
        set(TEST_RUN_DIR ${PROJECT_BINARY_DIR}/tests/${dir}/${_name})
    endif()
    file(MAKE_DIRECTORY ${TEST_RUN_DIR})

    set(INPUTFILE ${TEST_SRC_DIR}/input.dat)
    set(OUTFILE ${TEST_RUN_DIR}/output.dat)

    # Turn on psitest.pl if eligible
    # true/false have to be _lowercase_, CTest will otherwise spit out a BAD COMMAND error
    set(AUTOTEST false)
    if(labels)
        list(FIND labels "autotest" _index) 
        # If not found _index == -1; convert to FALSE
        if(PERL_FOUND AND (${_index} GREATER -1))
            set(AUTOTEST true)
        endif()
    endif()

    # AJ:
    # Sow/Reap tests need helpers to test the full routine
    # they should have the sowreap tag added 
    # this should extend to other Sow/Reap routines, I think 
    set(SOWREAP false)
    if(labels)
        list(FIND labels "sowreap" _index)
        if (${_index} GREATER -1)
            set(SOWREAP true)
        endif()
    endif()

    # Serial build
    add_test(NAME "${_name}"
      WORKING_DIRECTORY "${TEST_RUN_DIR}"
      COMMAND "${PYTHON_EXECUTABLE}" "${TESTEXE}" "${INPUTFILE}" "${LOGFILE}" "${AUTOTEST}" "${psi4_SOURCE_DIR}" "${SOWREAP}" "${OUTFILE}" "${PSIEXE}" "${PSIDATADIR}" "${PSILIB}"
    )

    if(labels)
        set_tests_properties(${_name} PROPERTIES LABELS "${labels}")
    endif()
endmacro()

# This macro is used to add a unit test using Google Unit Testing framework
macro(add_googletest test_name my_libraries external_libraries)
    get_filename_component(the_name ${test_name} NAME_WE)
    add_executable(${the_name}.x ${the_name}.cpp)
    add_dependencies(${the_name}.x googletest)
    target_link_libraries(${the_name}.x
                          ${my_libraries}
                          ${GTEST_LIBS_DIR}/libgtest.a
                          ${GTEST_LIBS_DIR}/libgtest_main.a
                          ${CMAKE_THREAD_LIBS_INIT}
                          ${external_libraries}
    )
    add_test(NAME ${the_name} COMMAND ${the_name}.x)
endmacro()

# This macro is used to copy a file containing reference values into the directory
# where the unit tests will be executed
macro(add_reference reference_file where)
    FILE(COPY ${reference_file} DESTINATION ${where}
         FILE_PERMISSIONS OWNER_READ OWNER_WRITE GROUP_READ
         WORLD_READ)
endmacro()

