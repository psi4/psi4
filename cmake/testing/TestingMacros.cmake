# This macro is used to add regression tests.
# Tests are identified by their name (one) and by labels (more than one)
# all tests with label short should run in less than 30 seconds sequential on a 2GHz CPU
macro(add_regression_test _name _labels)
    add_test(
        ${_name}
	WORKING_DIRECTORY ${TEST_RUN_DIR}
	COMMAND ${PYTHON_EXECUTABLE} ${TESTEXE} ${INPUTFILE} ${LOGFILE} ${AUTOTEST} ${PROJECT_SOURCE_DIR} ${OUTFILE} ${PSIEXE})
    if(NOT "${_labels}" STREQUAL "")
        set_tests_properties(${_name} PROPERTIES LABELS "${_labels}")
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

# This macro is used to add a unit test using Boost Unit Testing framework
macro(add_boosttest test_name)
    get_filename_component(the_name ${test_name} NAME_WE)
    add_executable(${the_name}.x ${the_name}.cpp)      	

    set(_my_libraries "${ARGV1}")
    set(_external_libraries "${ARGV2}")
    # Find threading library aka CMAKE_THREAD_LIBS_INIT
    find_package(Threads)

    # Building on more than one processor can result in race conditions,
    # since custom Boost can be built only on one processor!
    # We thus add this dependency to not get stuck.
    if(BUILD_CUSTOM_BOOST)
	    add_dependencies(${the_name}.x custom_boost)
    endif()	    
    target_link_libraries(${the_name}.x
                          ${_my_libraries}
			  ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
	    	          ${CMAKE_THREAD_LIBS_INIT}
                          ${_external_libraries}
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

