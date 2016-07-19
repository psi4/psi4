#Guard against in-source builds
macro(no_in_source)
if(${CMAKE_BINARY_DIR}==${CMAKE_SOURCE_DIR})
   message(WARNING "In source builds are prohibited.  Making a build directory")
   set(CMAKE_BINARY_DIR build)
endif()
endmacro(no_in_source)

#Wraps an option with a default other than ON/OFF
macro(option_with_default variable msge default)
if(NOT variable OR NOT ${variable})
   message(STATUS "No value for " ${variable} " specified using: " ${default})
   set(${variable} ${default} CACHE STRING ${msge})
endif()
endmacro(option_with_default)

#Common guts to adding a Psi4 library irrespective of bin vs. lib home
macro(general_add_library libname sources prefix dir)
   add_library(${libname} ${${sources}})
   set_target_properties(${libname} PROPERTIES 
       POSITION_INDEPENDENT_CODE ${BUILD_FPIC}
   )
   install(TARGETS ${libname} DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
   install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} 
      DESTINATION psi4/src/${dir}/${prefix}${libname}
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
macro(psi4_add_library libname sources)
   general_add_library(${libname} ${sources} lib lib ${ARGN}) 
endmacro(psi4_add_library libname sources)

#Adds a psi4 library that lives in bin
macro(psi4_add_binary libname sources)
   general_add_library(${libname} ${sources} "" bin ${ARGN})
endmacro(psi4_add_binary libname sources)

include(CheckCCompilerFlag)
include(CheckCXXCompilerFlag)

macro(add_C_flags FLAGS)
   unset(test_option)
   CHECK_C_COMPILER_FLAG("${FLAGS} -Werror" test_option)
   if(${test_option})
      set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${FLAGS}")
   endif()
endmacro()

macro(add_CXX_flags FLAGS)
   unset(test_option)
   CHECK_CXX_COMPILER_FLAG("${FLAGS} -Werror" test_option)
   if(${test_option})
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${FLAGS}")
   endif()
endmacro()

macro(add_flags FLAGS)
   add_C_flags(${FLAGS})
   add_CXX_flags(${FLAGS})
endmacro()

#Adds flags if compiler recognizes them and the option is true
macro(optional_add_flags option)
   set(Extra_Args ${ARGN})
   list(LENGTH Extra_Args nargs)
   if(${option} AND (${nargs} GREATER 0))
      list(GET Extra_Args 0 FLAGS)
      add_flags(${FLAGS})
   endif()
endmacro()

macro(option_with_flags option msg default)
    option(${option} ${msg} ${default})
    optional_add_flags(${option} ${ARGN})
endmacro()

#Macro so I don't have to look at a ton of if statements
macro(optional_plugin plugin_name)
string(TOUPPER ${plugin_name} PLUGIN_NAME)
if(${ENABLE_${PLUGIN_NAME}})
   find_package(${plugin_name} REQUIRED)
   set_property(GLOBAL APPEND PROPERTY PSI4_MODULES ${${PLUGIN_NAME}_LIBRARIES})
   add_definitions(-DENABLE_${PLUGIN_NAME})
endif()
endmacro(optional_plugin plugin_name)