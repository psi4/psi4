# Run the git version script, and update gitversion.h, if it's changed
add_custom_target(update_version
    COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/gitversion.py ${PROJECT_SOURCE_DIR}
    COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_CURRENT_BINARY_DIR}/gitversion.h.tmp ${CMAKE_CURRENT_BINARY_DIR}/gitversion.h
    COMMAND ${CMAKE_COMMAND} -E remove ${CMAKE_CURRENT_BINARY_DIR}/gitversion.h.tmp
    COMMENT "Generating Git info")

add_library(versioned_code OBJECT version.cc python.cc psi_start.cc)
add_dependencies(versioned_code update_version)
add_dependencies(versioned_code mints)

set(headers_list "")
# List of headers
list(APPEND headers_list script.h psi4.h)

# If you want to remove some headers specify them explictly here
if(DEVELOPMENT_CODE)
   list(REMOVE_ITEM headers_list "")
else()
   list(REMOVE_ITEM headers_list "")
endif()
# Sort alphabetically
list(SORT headers_list)

set(sources_list "")
# List of sources
list(APPEND sources_list export_psio.cc export_mints.cc psi_stop.cc export_functional.cc export_oeprop.cc export_plugins.cc export_blas_lapack.cc psi4.cc export_cubefile.cc export_benchmarks.cc export_efp.cc export_cubeprop.cc clean.cc create_new_plugin.cc script.cc set_memory.cc read_options.cc export_libparallel.cc)
list(APPEND sources_list $<TARGET_OBJECTS:versioned_code>)

# If you want to remove some sources specify them explictly here
if(DEVELOPMENT_CODE)
   list(REMOVE_ITEM sources_list "")
else()
   list(REMOVE_ITEM sources_list "")
endif()


get_property(from_src_bin GLOBAL PROPERTY PSILIB)
get_property(from_src_lib GLOBAL PROPERTY LIBLIST)
get_property(from_interfaces GLOBAL PROPERTY INTERFACES)
list(APPEND psi_libraries "${from_src_bin}" "${from_src_lib}" "${from_interfaces}")

# Prepare for linking executable
list(APPEND
     LINKLIBS "${LIBC_INTERJECT}"
     "${psi_libraries}" "${Boost_LIBRARIES}"
     "${LAPACK_LIBRARIES}" "${BLAS_LIBRARIES}"
     "${PYTHON_LIBRARY}" "${EXTERNAL_LIBS}"
     "${LIBUTIL_LIBRARIES}" "${LIBDL_LIBRARIES}"
     "${LIBRT_LIBRARIES}" "${LIBM_LIBRARIES}"
     )

# Executable psi4
add_executable(psi4 ${sources_list})
add_dependencies(psi4 update_version)
if(CUSTOM_BOOST_BUILD)
   add_dependencies(psi4 custom_boost)
endif()
target_link_libraries(psi4 ${LINKLIBS})
if(ENABLE_CHEMPS2)
    target_link_libraries(psi4 CHEMPS2::CHEMPS2)
endif()
if(ENABLE_PCMSOLVER)
  target_link_libraries(psi4 PCMSolver::PCMSolver)
endif()

# standalone python module psi4.so
add_executable(psi4so EXCLUDE_FROM_ALL ${sources_list})
if(CUSTOM_BOOST_BUILD)
   add_dependencies(psi4so custom_boost)
endif()
set_property(TARGET psi4so PROPERTY COMPILE_DEFINITIONS MAKE_PYTHON_MODULE)
set_property(TARGET psi4so PROPERTY LINK_FLAGS "-shared")
set_property(TARGET psi4so PROPERTY OUTPUT_NAME "psi4${CMAKE_EXECUTABLE_SUFFIX}.so")
target_link_libraries(psi4so "${LINKLIBS}")
