#
#  This macro sets up Google Test framework for unit testing
#
macro(setup_googletest)
	message(STATUS "Setting up Google Test Unit Testing Framework")
	configure_file(
    	${PROJECT_SOURCE_DIR}/src/utils/gtestPimpl.hpp.in
	${PROJECT_BINARY_DIR}/include/gtestPimpl.hpp)
	set(GTestCMakeArgs "${ExternalProjectCMakeArgs}")
	list(REMOVE_ITEM GTestCMakeArgs "-DBOOST_INCLUDEDIR=${Boost_INCLUDE_DIRS}" "-DBOOST_LIBRARYDIR=${Boost_LIBRARY_DIRS}")
        ExternalProject_Add(googletest
        	SVN_REPOSITORY http://googletest.googlecode.com/svn/trunk
               	DOWNLOAD_DIR ${PROJECT_SOURCE_DIR}
        	PREFIX ${PROJECT_SOURCE_DIR}/external
        	SOURCE_DIR ${PROJECT_SOURCE_DIR}/external/googletest
        	BINARY_DIR ${PROJECT_BINARY_DIR}/external/googletest-build
                STAMP_DIR ${PROJECT_BINARY_DIR}/external/googletest-stamp
        	TMP_DIR ${PROJECT_BINARY_DIR}/external/googletest-tmp
                INSTALL_DIR ${PROJECT_BINARY_DIR}/external
        	CMAKE_ARGS ${GTestCMakeArgs}
        	# Disable install step.
                INSTALL_COMMAND ""
        	)
        
        ExternalProject_Get_Property(googletest SOURCE_DIR)
        set(GTEST_INCLUDE_DIRS ${SOURCE_DIR}/include)
        ExternalProject_Get_Property(googletest BINARY_DIR)
        set(GTEST_LIBS_DIR ${BINARY_DIR})
        include_directories(${GTEST_INCLUDE_DIRS})
	find_package(Threads REQUIRED)
endmacro(setup_googletest)

#
#  This macro sets up Boost Test framework for unit testing
#  We are not doing anything special, as Boost Test doesn't require
#  anything more than including headers and linking to a library
#
macro(setup_boosttest)
	message(STATUS "Setting up Boost Test Unit Testing Framework")
endmacro(setup_boosttest)
