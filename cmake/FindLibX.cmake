#Code factorization for the Find various libraries for which CMake doesn't
#want to provide an already found version
#Syntax:
#
# find_lib_x(name mac_include linux_include [OtherName1 [ OtherName2 ...]])
#
#    name is w/o lib in front
#    mac_include is name of the MAC header file
#    linux_include is name of the Linux header file
#    Remaining arguments are used as alternative names for the library
#
macro(find_lib_x lib_name mac_include linux_include)
   string(TOUPPER ${lib_name} LIB_NAME)
   if(LIB${LIB_NAME}_INCLUDE_DIR)
      # Already in cache, be silent
      set(LIB${LIB_NAME}_FIND_QUIETLY TRUE)
   endif()
   if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
      find_path(LIB${LIB_NAME}_INCLUDE_DIR ${mac_include}
   	  PATHS
   	  /usr/include
   	  /usr/local/include
      )
   else()
      # Should probably be generalized for Windows
      find_path(LIB${LIB_NAME}_INCLUDE_DIR ${linux_include}
   	  PATHS
   	  /usr/include
   	  /usr/local/include
      )
   endif()
   
   list(APPEND LIB${LIB_NAME}_NAMES ${lib_name} lib${lib_name})
   set(Extra_Args "${ARGN}")
   foreach(name_i IN LISTS Extra_Args)
      list(APPEND LIB${LIB_NAME}_NAMES ${other_name} lib${other_name})
   endforeach()
   find_library(LIB${LIB_NAME}_LIBRARY 
	     NAMES ${LIB${LIB_NAME}_NAMES}
	     PATHS 
	     /usr/lib
	     /usr/local/lib
	     /lib
   )
   include(FindPackageHandleStandardArgs)
   find_package_handle_standard_args(lib${lib_name} DEFAULT_MSG 
                                  LIB${LIB_NAME}_LIBRARY 
                                  LIB${LIB_NAME}_INCLUDE_DIR
   )
   if(LIB${LIB_NAME}_FOUND)
      set(LIB${LIB_NAME}_LIBRARIES ${LIB${LIB_NAME}_LIBRARY})
   else(LIB${LIB_NAME}_FOUND)
      set(LIB${LIB_NAME}_LIBRARIES)
   endif(LIB${LIB_NAME}_FOUND)

   mark_as_advanced(LIB${LIB_NAME}_LIBRARY LIB${LIB_NAME}_INCLUDE_DIR)
endmacro()
