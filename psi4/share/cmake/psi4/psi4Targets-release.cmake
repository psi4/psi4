#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "psi4::core" for configuration "Release"
set_property(TARGET psi4::core APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(psi4::core PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/psi4/core.cpython-312-x86_64-linux-gnu.so"
  IMPORTED_SONAME_RELEASE "core.cpython-312-x86_64-linux-gnu.so"
  )

list(APPEND _cmake_import_check_targets psi4::core )
list(APPEND _cmake_import_check_files_for_psi4::core "${_IMPORT_PREFIX}/lib/psi4/core.cpython-312-x86_64-linux-gnu.so" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
