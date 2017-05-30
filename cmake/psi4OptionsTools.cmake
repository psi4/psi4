###This file contains functions used throughout the Psi4 build.  Like source
###code, the build system should be factored and common code extracted out into
###functions/macros.  If you find repetitive code throughout the build scripts
###this is the place to add it (make sure you document it too).

#Macro for printing an option in a consistent manner
#
#Syntax: print_option(<option to print> <was specified>)
#
macro(print_option variable default)
if(NOT DEFINED ${variable} OR "${${variable}}" STREQUAL "")
message(STATUS "Setting (unspecified) option ${variable}: ${default}")
else()
message(STATUS "Setting option ${variable}: ${${variable}}")
endif()
endmacro()

# Wraps an option with default ON/OFF. Adds nice messaging to option()
#
#Syntax: option_with_print(<option name> <description> <default value>)
#
macro(option_with_print variable msge default)
print_option(${variable} ${default})
option(${variable} ${msge} ${default})
endmacro(option_with_print)

#Wraps an option with a default other than ON/OFF and prints it
#NOTE: Can't combine with above b/c CMake handles ON/OFF options specially
#NOTE2: CMAKE_BUILD_TYPE (and other CMake variables) are always defined so need
#       to further check for if they are the NULL string.  This is also why we
#       need the force
#
#Syntax: option_with_default(<option name> <description> <default value>)
#
macro(option_with_default variable msge default)
print_option(${variable} ${default})
if(NOT DEFINED ${variable} OR "${${variable}}" STREQUAL "")
set(${variable} ${default} CACHE STRING ${msge} FORCE)
endif()
endmacro(option_with_default)

# Common guts to adding a Psi4 library irrespective of bin vs. lib home
#
# Syntax psi4_add_module(<lib or bin> <library name> <CMake list of sources> <dependencies>)
#
macro(psi4_add_module binlib libname sources)

    set(current_sources ${${sources}};)
    list(SORT current_sources)

    add_library(${libname} STATIC ${current_sources})
    set_target_properties(${libname} PROPERTIES POSITION_INDEPENDENT_CODE ${BUILD_FPIC})

    # library modules get their headers installed
    if((${binlib} MATCHES lib) OR (${binlib} MATCHES binlib))
        install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}/psi4
                FILES_MATCHING PATTERN "*.h" PATTERN "*.hpp" PATTERN "*.i")
    endif()

    # binary modules explicitly compiled into psi4.so
    if((${binlib} MATCHES bin) OR (${binlib} MATCHES binlib))
        set_property(GLOBAL APPEND PROPERTY BINLIST ${libname})
    endif()

    set(depend_name "${ARGN}")
    foreach(name_i IN LISTS depend_name)
        target_link_libraries(${libname} PRIVATE ${name_i})
    endforeach()
    target_link_libraries(${libname} PRIVATE pybind11::module)
    target_link_libraries(${libname} PRIVATE tgt::lapack)
endmacro()

include(CheckCCompilerFlag)
include(CheckCXXCompilerFlag)
if(CMAKE_Fortran_COMPILER)
    include(CheckFortranCompilerFlag)  # CMake >= 3.3, so local copy in cmake/
endif()

#The guts of the next two functions, use the wrappers please
#
#Syntax: add_C_or_CXX_flags(<True for C, False for CXX>)
#
# Note: resist adding -Werror to the check_X_compiler_flag calls,
#   as (i) the flag for Intel is actually -diag-error warn, (ii)
#   Intel ifort doesn't define -Werror, and (iii) passing it
#   changes REQUIRED_DEFINITIONS.
macro(add_C_or_CXX_flags is_C)
set(CMAKE_REQUIRED_QUIET_SAVE ${CMAKE_REQUIRED_QUIET})
   set(CMAKE_REQUIRED_QUIET ON)
   set(flags_to_try "${ARGN}")
   foreach(flag_i IN LISTS flags_to_try ITEMS -brillig)
      if(${flag_i} STREQUAL "-brillig")
         message(WARNING "Option unfulfilled as none of ${flags_to_try} valid")
         break()
      endif()
      unset(test_option CACHE)
      if(${is_C} EQUAL 0)
          CHECK_C_COMPILER_FLAG("${flag_i}" test_option)
          set(description_to_print CMAKE_C_FLAGS)
      elseif(${is_C} EQUAL 1)
          CHECK_CXX_COMPILER_FLAG("${flag_i}" test_option)
          set(description_to_print CMAKE_CXX_FLAGS)
      elseif(${is_C} EQUAL 2)
          CHECK_Fortran_COMPILER_FLAG("${flag_i}" test_option)
          set(description_to_print CMAKE_Fortran_FLAGS)
      endif()
      set(msg_base "Performing Test ${description_to_print} [${flag_i}] -")
      if(${test_option})
        set(${description_to_print} "${${description_to_print}} ${flag_i}")
        if(NOT CMAKE_REQUIRED_QUIET_SAVE)
           message(STATUS  "${msg_base} Success, Appending")
        endif()
        break()
      else()
        if(NOT CMAKE_REQUIRED_QUIET_SAVE)
           message(STATUS "${msg_base} Failed")
        endif()
      endif()
   endforeach()
   set(CMAKE_REQUIRED_QUIET ${CMAKE_REQUIRED_QUIET_SAVE})
endmacro()



#Checks if C flags are valid, if so adds them to CMAKE_C_FLAGS
#Input should be a list of flags to try.  If two flags are to be tried together
#enclose them in quotes, e.g. "-L/path/to/dir -lmylib" is tried as a single
#flag, whereas "-L/path/to/dir" "-lmylib" is tried as two separate flags.
#The first list item to succeed is added to CMAKE_C_FLAGS, then try loop
#breaks. Warning issued if no flags in list succeed.
#
#
#Syntax: add_C_flags(<flags to add>)
#
macro(add_C_flags)
   add_C_or_CXX_flags(0 ${ARGN})
endmacro()

#Checks if CXX flags are valid, if so adds them to CMAKE_CXX_FLAGS
#See add_C_flags for more info on syntax
#
#Syntax: add_CXX_flags(<flags to add>)
#
macro(add_CXX_flags)
    add_C_or_CXX_flags(1 ${ARGN})
endmacro()

#Checks if Fortran flags are valid, if so adds them to CMAKE_Fortran_FLAGS
#See add_C_flags for more info on syntax
#
#Syntax: add_Fortran_flags(<flags to add>)
#
macro(add_Fortran_flags)
    add_C_or_CXX_flags(2 ${ARGN})
endmacro()

#Macro for adding flags common to both C and CXX, if the compiler supports them
#
#Syntax: add_flags(<flags to add>)
#
macro(add_flags FLAGS)
    get_property(languages GLOBAL PROPERTY ENABLED_LANGUAGES)
    list(FIND languages "C" _index_c)
    list(FIND languages "CXX" _index_cxx)
    list(FIND languages "Fortran" _index_fortran)
    if (${_index_c} GREATER -1)
        add_C_flags(${FLAGS})
    endif()
    if (${_index_cxx} GREATER -1)
        add_CXX_flags(${FLAGS})
    endif()
    if (${_index_fortran} GREATER -1)
        add_Fortran_flags(${FLAGS})
    endif()
endmacro()

#Defines an option that if enabled turns on some compiler flags
#
#Syntax: option_with_flags(<option> <description> <default value> <flags>)
#
macro(option_with_flags option msg default)
    print_option(${option} ${default})
    option(${option} ${msg} ${default})
    if(${${option}})
       add_flags("${ARGN}")
    endif()
endmacro()

#Macro so I don't have to look at a ton of if statements for adding each plugin
#
#Syntax: optional_plugin(<plugin name>)
#
macro(optional_plugin plugin_name)
string(TOUPPER ${plugin_name} PLUGIN_NAME)
if(${ENABLE_${PLUGIN_NAME}})
   find_package(${plugin_name} REQUIRED)
   set_property(GLOBAL APPEND PROPERTY PSI4_MODULES ${${PLUGIN_NAME}_LIBRARIES})
   add_definitions(-DENABLE_${PLUGIN_NAME})
else()
   add_library(${plugin_name} INTERFACE)
endif()
endmacro(optional_plugin plugin_name test_names)

#Macro for adding a skeleton plugin
macro(add_skeleton_plugin PLUG TEMPLATE TESTLABELS)
    set(CCSD "${CMAKE_CURRENT_SOURCE_DIR}")
    set(CCBD "${CMAKE_CURRENT_BINARY_DIR}")
    set(PSIEXE ${STAGED_INSTALL_PREFIX}/${CMAKE_INSTALL_BINDIR}/psi4)
    set(DIR_2_PASS ${CMAKE_PREFIX_PATH} ${STAGED_INSTALL_PREFIX})

    add_custom_target(plugin_${PLUG}
        ALL
        DEPENDS psi4-core
        COMMAND ${CMAKE_COMMAND} -E remove_directory ${CCBD}/${PLUG}
        COMMAND ${PSIEXE} --plugin-name ${PLUG} --plugin-template ${TEMPLATE}
        COMMAND ${CMAKE_COMMAND} -E chdir "${CCBD}/${PLUG}" cmake -C ${STAGED_INSTALL_PREFIX}/share/cmake/psi4/psi4PluginCache.cmake "-DCMAKE_PREFIX_PATH=${DIR_2_PASS}" .
        COMMAND ${CMAKE_COMMAND} -E chdir "${CCBD}/${PLUG}" ${CMAKE_MAKE_PROGRAM}
        COMMAND ${CMAKE_COMMAND} -E create_symlink ${CCBD}/${PLUG}/input.dat ${CCSD}/input.dat
        COMMAND ${CMAKE_COMMAND} -E create_symlink "${PLUG}/${PLUG}.so" "${PLUG}.so"
        COMMAND ${CMAKE_COMMAND} -E create_symlink "${PLUG}/__init__.py" "__init__.py"
        COMMAND ${CMAKE_COMMAND} -E create_symlink "${PLUG}/pymodule.py" "pymodule.py"
        COMMENT "Build ${PLUG} example plugin"
        VERBATIM)

    include(TestingMacros)
    add_regression_test(${PLUG} "${TESTLABELS}")
endmacro()


