#
# CTest script for automated project testing and submission to CDash
#
# Written by Jonas Juselius <jonas.juselius@uit.no>, November 2010
#
# This script will run a series of tests and submit them to the project 
# dashboard (configured in the project directory). By default it runs
# test for the Experimental model, on the master branch with a Debug build.
#
# The script accepts arguments in the following manner:
#
# $ ctest -C Debug -S CDashTestBranch.cmake,branch=master,model=Nightly,...
#
# Valid options are:
#   branch=name
#   model=(Experimental|Nightly|Continuous)
#   site=name
#   toolchain_name=name (e.g. Intel, GNU)
#   cc=C compiler
#   cxx=C++ compiler
#   fc=Fortran compiler
#   mpi=(on|off)
#   openmp=(on|off)
#   memcheck=(on|off)
#   coverage(on|off)
#
# The actual tests runs are defined at the end of the script
#

# -------------------------------------------------------------------------
# -- User configuration
# -------------------------------------------------------------------------

if(DEFINED ENV{PROJECT_NAME})
    set (PROJECT_NAME $ENV{PROJECT_NAME})
endif()

if(NOT DEFINED PROJECT_NAME)
    message(FATAL_ERROR "PROJECT_NAME not set!")
endif()

if(DEFINED ENV{PROJECT_REPOSITORY})
    set (PROJECT_REPOSITORY $ENV{PROJECT_REPOSITORY})
endif()

if(NOT DEFINED PROJECT_REPOSITORY)
    message(FATAL_ERROR "PROJECT_REPOSITORY not set!")
endif()

if(DEFINED ENV{DASHBOARD_DIR})
    set (DASHBOARD_DIR $ENV{DASHBOARD_DIR}/${PROJECT_NAME})
else()
    set (DASHBOARD_DIR "/tmp/${PROJECT_NAME}")
endif()

set(CTEST_COMMAND ctest)
set(CTEST_CMAKE_GENERATOR "Unix Makefiles")

# -------------------------------------------------------------------------
# -- Initialize the CTest environment
# -------------------------------------------------------------------------
if(DEFINED ENV{CMAKE_CONFIG_TYPE})
    set(BUILD_TYPE $ENV{CMAKE_CONFIG_TYPE})
else()
    set(BUILD_TYPE Debug)
endif()

find_program(CTEST_GIT_COMMAND NAMES git)
find_program(CTEST_COVERAGE_COMMAND NAMES gcov)
find_program(CTEST_MEMORYCHECK_COMMAND NAMES valgrind)
find_program(UNAME NAMES uname)
find_program(HOSTNAME_COMMAND NAMES hostname)

set(CTEST_UPDATE_COMMAND ${CTEST_GIT_COMMAND})

include(ProcessorCount)
ProcessorCount(NCORES)
if(NOT NCORES EQUAL 0)
    set(CTEST_BUILD_FLAGS -j${NCORES})
    #set(ctest_test_args ${ctest_test_args} PARALLEL_LEVEL ${N})
else()
    set(NCORES 1)
endif()

# --- Defaults
set(BRANCH master)
set(TEST_MODEL Experimental)
set(ENABLE_MPI FALSE)
set(ENABLE_OMP FALSE)
set(ENABLE_MEMCHECK TRUE)
set(ENABLE_COVERAGE TRUE)
execute_process(COMMAND ${HOSTNAME_COMMAND} 
    OUTPUT_VARIABLE CTEST_SITE
    ERROR_QUIET
    )
string(STRIP ${CTEST_SITE} CTEST_SITE)
# -------------------------------------------------------------------------
# -- Given a variable VAR and a string FOO=bar -> set (VAR "bar")
# -------------------------------------------------------------------------
macro(set_var_from_str var str)
string(REPLACE "=" " " arglist ${str})
separate_arguments(arglist)
list(GET arglist 1 val)
set(${var} "${val}")
unset(arglist)
unset(val)
endmacro()

# -------------------------------------------------------------------------
# -- Parse the command line args
# -------------------------------------------------------------------------
macro(parse_args)
string(REPLACE "," " " arglist ${CTEST_SCRIPT_ARG})
separate_arguments(arglist)
foreach(arg ${arglist})
    if("${arg}" MATCHES "branch=")
        set_var_from_str(BRANCH ${arg})
    endif()

    if("${arg}" MATCHES "model=")
        set_var_from_str(TEST_MODEL ${arg})
    endif()

    if("${arg}" MATCHES "toolchain_name=")
        set_var_from_str(TOOLCHAIN_NAME ${arg})
    endif()

    if("${arg}" MATCHES "cxx=")
        set_var_from_str(CXX ${arg})
    endif()

    if("${arg}" MATCHES "cc=")
        set_var_from_str(CC ${arg})
    endif()

    if("${arg}" MATCHES "fc=")
        set_var_from_str(FC ${arg})
    endif()

    if("${arg}" MATCHES "mpi=")
        set_var_from_str(ENABLE_MPI ${arg})
    endif()

    if("${arg}" MATCHES "openmp=")
        set_var_from_str(ENABLE_OMP ${arg})
    endif()

    if("${arg}" MATCHES "memcheck=(no|off|false|False)")
        set (ENABLE_MEMCHECK FALSE)
    endif()

    if("${arg}" MATCHES "coverage=(no|off|false|False)")
        set (ENABLE_COVERAGE FALSE)
    endif()

    if("${arg}" MATCHES "site=")
        set_var_from_str(CTEST_SITE ${arg})
    endif()
endforeach()
unset(arglist)
unset(arg)
endmacro()

parse_args()

# 
# -- Initialize directories
# 
set(CTEST_SOURCE_DIRECTORY 
    ${DASHBOARD_DIR}/source/${BRANCH}
    )
set(CTEST_BINARY_DIRECTORY 
    ${DASHBOARD_DIR}/build/${BRANCH}
    )
    endif()

if(NOT EXISTS ${CTEST_SOURCE_DIRECTORY})
    execute_process(COMMAND ${CTEST_GIT_COMMAND} 
        clone ${PROJECT_REPOSITORY} ${CTEST_SOURCE_DIRECTORY}
        OUTPUT_QUIET
        )
    if (NOT (${BRANCH} STREQUAL master))
        execute_process(COMMAND ${CTEST_GIT_COMMAND} checkout 
            -b ${BRANCH} origin/${BRANCH}
            WORKING_DIRECTORY ${CTEST_SOURCE_DIRECTORY}
            OUTPUT_QUIET
            )
    endif()
endif()

# -------------------------------------------------------------------------
# -- Configure the build name based on the OS and CPU + additional args
# -------------------------------------------------------------------------
macro(get_build_name name)
execute_process(COMMAND ${UNAME} -s OUTPUT_VARIABLE osname)
string(STRIP ${osname} osname)
execute_process(COMMAND ${UNAME} -m OUTPUT_VARIABLE cpu)
string(STRIP ${cpu} cpu)

set(${name} "${osname}-${cpu}")
foreach(arg ${ARGN})
    set(${name} "${${name}}-${arg}")
endforeach()
unset(osname)
unset(cpu)
unset(arg)
endmacro()

# -------------------------------------------------------------------------
# -- Append to a variable
# -------------------------------------------------------------------------
macro(append_to_var name)
foreach(arg ${ARGN})
    set(${name} "${${name}}${arg}")
endforeach()
unset(arg)
endmacro()

# -------------------------------------------------------------------------
# -- Given a string FOO=bar, export it to the environment
# -------------------------------------------------------------------------
function(set_env_from_string str)
string(REPLACE "=" " " envlist ${str})
separate_arguments(envlist)
list(GET envlist 0 envvar)
list(GET envlist 1 envval)
set(ENV{${envvar}} ${envval})
unset(envlist)
unset(envvar)
endfunction()

# -------------------------------------------------------------------------
# -- Run the shell script 'setup' which produces the relevant CMake flags
# -- and environment variables for the build. Upon exit environment variables
# -- have been set, and CTEST_COVERAGE_COMMAND is defined.
# -------------------------------------------------------------------------
macro(setup_cmake_env)
get_build_name(CTEST_BUILD_NAME ${BRANCH})

set (SETUP_FLAGS --enable-tests --show)

if (${BUILD_TYPE} STREQUAL Debug)
    set (SETUP_FLAGS ${SETUP_FLAGS} --debug)
    append_to_var(CTEST_BUILD_NAME "-dbg")
elseif (${BUILD_TYPE} STREQUAL Release)
    set (SETUP_FLAGS ${SETUP_FLAGS} --release)
    append_to_var(CTEST_BUILD_NAME "-rel")
endif()
        
if (ENABLE_MPI)
    set (SETUP_FLAGS ${SETUP_FLAGS} --enable-mpi)
    append_to_var(CTEST_BUILD_NAME "-mpi")
endif()

if (ENABLE_OMP)
    set (SETUP_FLAGS ${SETUP_FLAGS} --enable-openmp)
    append_to_var(CTEST_BUILD_NAME "-omp")
endif()

if (ENABLE_COVERAGE)
    set (SETUP_FLAGS ${SETUP_FLAGS} --coverage)
endif()

if (DEFINED TOOLCHAIN_NAME)
    append_to_var(CTEST_BUILD_NAME "-${TOOLCHAIN_NAME}")
endif()

if (DEFINED CTEST_SITE)
    set (SETUP_FLAGS ${SETUP_FLAGS} --host=${CTEST_SITE})
endif()

if (DEFINED CC)
    set (SETUP_FLAGS ${SETUP_FLAGS} --cc=${CC})
endif()

if (DEFINED CXX)
    set (SETUP_FLAGS ${SETUP_FLAGS} --cxx=${CXX})
endif()

if (DEFINED FC)
    set (SETUP_FLAGS ${SETUP_FLAGS} --fc=${FC})
endif()

execute_process(
    COMMAND ${CTEST_SOURCE_DIRECTORY}/setup ${SETUP_FLAGS}
    OUTPUT_VARIABLE setup_list
    )
separate_arguments(setup_list)
foreach(arg ${setup_list})
    if("${arg}" MATCHES "-.+=.*")
        set(CTEST_CONFIGURE_COMMAND "${CTEST_CONFIGURE_COMMAND} ${arg}")
    elseif("${arg}" MATCHES ".+=.*")
        set_env_from_string(${arg})
    else()
        set(CTEST_CONFIGURE_COMMAND "${CTEST_CONFIGURE_COMMAND} ${arg}")
    endif()
endforeach()

set (CTEST_CONFIGURE_COMMAND 
    "${CTEST_CONFIGURE_COMMAND} ${CTEST_SOURCE_DIRECTORY}"
    )
unset(arg)
unset(setup_list)
endmacro()

# -------------------------------------------------------------------------
# -- Run the dashboard
# -------------------------------------------------------------------------
macro(run_dashboard)
message("${CTEST_SOURCE_DIRECTORY}/setup ${SETUP_FLAGS}")
ctest_configure()
ctest_build()
ctest_test(PARALLEL_LEVEL ${NCORES})
if(ENABLE_COVERAGE AND CTEST_COVERAGE_COMMAND)
    ctest_coverage()
endif ()
if(ENABLE_MEMCHECK AND CTEST_MEMORYCHECK_COMMAND)
    ctest_memcheck(PARALLEL_LEVEL ${NCORES})
endif ()
endmacro()


# -------------------------------------------------------------------------
# -- Run the tests
# -------------------------------------------------------------------------

setup_cmake_env()
ctest_read_custom_files(${CTEST_BINARY_DIRECTORY})
ctest_empty_binary_directory(${CTEST_BINARY_DIRECTORY})
ctest_start(${TEST_MODEL})
ctest_update()
run_dashboard()
ctest_submit()

# vim:et:sw=4:ts=4:
