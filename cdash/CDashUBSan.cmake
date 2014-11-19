#
# CTest script for undefined behaviour sanitizer (UBSan) memcheck build and submission to dashboard 
#
# Written by Roberto Di Remigio November 2014 
#

set(PROJECT_NAME "Psi")
set(PROJECT_REPOSITORY "git@github.com:robertodr/psi4.git") 
set(CTEST_COMMAND ctest)
set(CTEST_CMAKE_GENERATOR "Unix Makefiles")
set(CTEST_BUILD_NAME "RDR-gcc4.9.1-ubsan")
set(CTEST_SITE "stallo.uit.no")

find_program(CTEST_GIT_COMMAND NAMES git)
set(CTEST_UPDATE_COMMAND ${CTEST_GIT_COMMAND})


include(ProcessorCount)
ProcessorCount(NCORES)
if(NOT NCORES EQUAL 0)
    set(CTEST_BUILD_FLAGS -j${NCORES})
else()
    set(NCORES 1)
endif()

# Defaults
set(BRANCH master)
set(TEST_MODEL Experimental)

# Initialize directories
set(CTEST_SOURCE_DIRECTORY $ENV{CTEST_SOURCE_DIRECTORY})
set(CTEST_BINARY_DIRECTORY $ENV{CTEST_BINARY_DIRECTORY})
# Fresh clone of repository and checkout right branch
if(NOT EXISTS ${CTEST_SOURCE_DIRECTORY})
    execute_process(COMMAND ${CTEST_GIT_COMMAND} 
        clone ${PROJECT_REPOSITORY} ${CTEST_SOURCE_DIRECTORY}
        OUTPUT_QUIET
        )
    if(NOT (${BRANCH} STREQUAL master))
        execute_process(COMMAND ${CTEST_GIT_COMMAND} checkout 
            -b ${BRANCH} origin/${BRANCH}
            WORKING_DIRECTORY ${CTEST_SOURCE_DIRECTORY}
            OUTPUT_QUIET
            )
    endif()
endif()

# Run the setup script
macro(setup_cmake_env)
    set(CTEST_CONFIGURE_COMMAND "cmake  -DCMAKE_Fortran_COMPILER=gfortran -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ -DENABLE_MPI=OFF -DENABLE_SGI_MPT=OFF -DENABLE_OMP=ON -DENABLE_VECTORIZATION=OFF -DENABLE_CSR=OFF -DENABLE_SCALAPACK=OFF -DENABLE_SCALASCA=OFF -DENABLE_UNIT_TESTS=OFF -DENABLE_STATIC_LINKING=OFF -DENABLE_PLUGINS=ON -DENABLE_LIBERD=OFF -DENABLE_GPU_DFCC=OFF -DENABLE_DUMMY_PLUGIN=OFF -DENABLE_CXX11_SUPPORT=ON -DLIBINT_OPT_AM=5 -DENABLE_UBSAN=ON -DCMAKE_INSTALL_PREFIX=/usr/local/psi4 -DCMAKE_BUILD_TYPE=debug ${CTEST_SOURCE_DIRECTORY}")
endmacro()    

# Run the dashboard
macro(run_dashboard)
    ctest_configure()
    ctest_build()
    ctest_test(PARALLEL_LEVEL ${NCORES})
    ctest_memcheck(PARALLEL_LEVEL ${NCORES})
endmacro()

# This is the actual run
setup_cmake_env()
ctest_read_custom_files(${CTEST_BINARY_DIRECTORY})
ctest_empty_binary_directory(${CTEST_BINARY_DIRECTORY})
set(CTEST_MEMORYCHECK_TYPE "UndefinedBehaviorSanitizer")
ctest_start(${TEST_MODEL})
ctest_update()
run_dashboard()
ctest_submit()

# vim:et:sw=4:ts=4:
