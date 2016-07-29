###This file contains functions used throughout the Psi4 build.  Like source
###code, the build system should be factored and common code extracted out into
###functions/macros.  If you find repetitive code throughout the build scripts
###this is the place to add it (make sure you document it too).

#Guard against in-source builds
#
# Syntax no_in_source()
#
macro(no_in_source)
if(${CMAKE_BINARY_DIR}==${CMAKE_SOURCE_DIR})
   message(WARNING "In-source builds are prohibited. Making a build directory")
   set(CMAKE_BINARY_DIR build)
   file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR})
endif()
endmacro(no_in_source)

#Macro for printing an option in a consistent manner
#
#Syntax: print_option(<option to print> <was specified>)
#
macro(print_option variable default)
if(NOT DEFINED ${variable})
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
#
#Syntax: option_with_default(<option name> <description> <default value>)
#
macro(option_with_default variable msge default)
print_option(${variable} ${default})
if(NOT DEFINED ${variable})
   set(${variable} ${default} CACHE STRING ${msge})
endif()
endmacro(option_with_default)


#Common guts to adding a Psi4 library irrespective of bin vs. lib home
#NOTE: list of sources is a CMake list
#
#Syntax general_add_library(<library name>, <list of sources>, <lib or bin>,
#                           <dependencies>)
#
macro(general_add_library libname sources dir)
   #TODO: Switch to OBJECT library?  Simplifies this macro...
   if(${dir} MATCHES lib)
      set(prefix lib)
   endif()
   add_library(${libname} ${${sources}})
   set_target_properties(${libname} PROPERTIES 
       POSITION_INDEPENDENT_CODE ${BUILD_FPIC}
   )
   install(TARGETS ${libname} DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR})
   install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} 
      DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}/psi4/src/${dir}/${prefix}${libname}
      FILES_MATCHING PATTERN "*.h" PATTERN "*.hpp")
   set_property(GLOBAL APPEND PROPERTY LIBLIST ${libname})
   set(depend_name "${ARGN}")
   foreach(name_i IN LISTS depend_name)
      target_link_libraries(${libname} INTERFACE ${name_i})
   endforeach()
   target_include_directories(${libname} PUBLIC ${Boost_INCLUDE_DIRS} 
                                                ${LIBDERIV_INCLUDE_DIRS})
endmacro(general_add_library libname sources prefix dir)

#Adds a psi4 library that lives in lib
#
#Syntax: psi4_add_library(<library name> <list of sources> <dependencies>)
#
macro(psi4_add_library libname sources)
   general_add_library(${libname} ${sources} lib ${ARGN}) 
endmacro(psi4_add_library libname sources)

#Adds a psi4 library that lives in bin
#
#Syntax: psi4_add_binary(<library name> <list of sources> <dependencies>
#
macro(psi4_add_binary libname sources)
   general_add_library(${libname} ${sources} bin ${ARGN})
endmacro(psi4_add_binary libname sources)

include(CheckCCompilerFlag)
include(CheckCXXCompilerFlag)

#Checks if C flags are valid, if so adds them to CMAKE_C_FLAGS
#Input should be a list of flags to try.  If two flags are to be tried together
#enclose them in quotes, e.g. "-L/path/to/dir -lmylib" is tried as a single
#flag, whereas "-L/path/to/dir" "-lmylib" is tried as two separate flags
#
#Syntax: add_C_flags(<flags to add>)
#
macro(add_C_flags)
   set(flags_to_try "${ARGN}")
   foreach(flag_i IN LISTS flags_to_try)
      unset(test_option)
      CHECK_C_COMPILER_FLAG("${flag_i} -Werror" test_option)
      if(${test_option})
         set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${FLAGS}")
         message(STATUS "Appending CMAKE_C_FLAGS: ${FLAGS}")
      endif()
   endforeach()
endmacro()

#Checks if CXX flags are valid, if so adds them to CMAKE_CXX_FLAGS
#See add_C_flags for more info on syntax
#
#Syntax: add_CXX_flags(<flags to add>)
#
macro(add_CXX_flags)
   set(flags_to_try "${ARGN}")
   foreach(flag_i IN LISTS flags_to_try)
      unset(test_option)
      CHECK_CXX_COMPILER_FLAG("${FLAGS} -Werror" test_option)
      if(${test_option})
         set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${FLAGS}")
         message(STATUS "Appending CMAKE_CXX_FLAGS: ${FLAGS}")
      endif()
   endforeach()
endmacro()

#Macro for adding flags common to both C and CXX, if the compiler supports them
#
#Syntax: add_flags(<flags to add>)
#
macro(add_flags FLAGS)
   add_C_flags(${FLAGS})
   add_CXX_flags(${FLAGS})
endmacro()

#Adds flags if compiler recognizes them and the option is true.  Flags are not
#necessarily governed by an option, e.g. restrict keyword
#
#Syntax: optional_add_flags(<option> <flags to add>)
#
macro(optional_add_flags option)
   set(Extra_Args ${ARGN})
   list(LENGTH Extra_Args nargs)
   if(${option} AND (${nargs} GREATER 0))
      list(GET Extra_Args 0 FLAGS)
      add_flags(${FLAGS})
   endif()
endmacro()

#Defines an option that if enabled turns on some compiler flags
#
#Syntax: option_with_flags(<option> <description> <default value> <flags>)
#
macro(option_with_flags option msg default)
    print_option(${option} ${default})
    option(${option} ${msg} ${default})
    optional_add_flags(${option} ${ARGN})
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
endmacro(optional_plugin plugin_name)