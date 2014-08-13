#This fxn just generates a list of directories in the current directory
MACRO(SUBDIRLIST result curdir)
  FILE(GLOB children RELATIVE ${curdir} ${curdir}/*)
  SET(dirlist "")
  FOREACH(child ${children})
    IF(IS_DIRECTORY ${curdir}/${child})
        LIST(APPEND dirlist "${child}")
    ENDIF()
  ENDFOREACH()
  SET(${result} ${dirlist})
ENDMACRO()

SUBDIRLIST(SUBDIRS ${CMAKE_CURRENT_SOURCE_DIR})
 
#RMR---It's important that only real plugins have CMakeLists.txt because that
# is how this next function determines what plugins to link against t
FOREACH(subdir ${SUBDIRS})
  set(file2look4 "${CMAKE_CURRENT_SOURCE_DIR}/${subdir}/CMakeLists.txt")
  message("Checking ${file2look4} for CMakeLists.txt")
  if(EXISTS ${file2look4})
    set(PLUGINS ${PLUGINS} ${subdir} PARENT_SCOPE)
    add_subdirectory(${subdir})
  endif()
ENDFOREACH()
