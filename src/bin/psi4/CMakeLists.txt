include(GlobUtils)

# Run the git version script, and update gitversion.h, if it's changed
add_custom_command(
    COMMAND ${PYTHON_EXECUTABLE} gitversion.py
    COMMAND ${CMAKE_COMMAND} -E copy_if_different gitversion.h.tmp gitversion.h
    OUTPUT gitversion.h
    COMMENT "Generating Git info"
)

add_custom_command(
    COMMAND ${CMAKE_COMMAND} -E copy_if_different ${PROJECT_SOURCE_DIR}/src/bin/psi4/version.cc
                                                  ${PROJECT_BINARY_DIR}/src/bin/psi4/version.cc
    OUTPUT version.cc
    DEPENDS gitversion.h
)

set(headers_list "")
# Get headers with .h extension
get_h_headers(headers_list)
# If you want to remove some headers specify them explictly here
if(DEVELOPMENT_CODE)
   list(REMOVE_ITEM headers_list "")
else()
   list(REMOVE_ITEM headers_list "")
endif()
# Sort alphabetically
list(SORT headers_list)

set(sources_list "")
# Get sources with .cc extension
get_cc_sources(sources_list)
# If you want to remove some sources specify them explictly here
if(DEVELOPMENT_CODE)
   list(REMOVE_ITEM sources_list "")
else()
   list(REMOVE_ITEM sources_list "")
endif()

# Write list of files to be passed to cloc for counting lines of code.
# Only files that are actually compiled are counted.
set(to_count "${sources_list}" "${headers_list}")
write_to_cloc_list("${to_count}")

# Prepare for linking executable
list(APPEND LINKLIBS "${PSILIB}" "${Boost_LIBRARIES}" "${LAPACK_LIBRARIES}" "${BLAS_LIBRARIES}" "${PYTHON_LIBRARIES}" dl)

# Executable psi4
add_executable(psi4 ${sources_list})
if(CUSTOM_BOOST_BUILD)
   add_dependencies(psi4 custom_boost)
endif()
target_link_libraries(psi4 ${LINKLIBS})

# standalone python module psi4.so
add_executable(psi4so EXCLUDE_FROM_ALL ${sources_list})
set_property(TARGET psi4so PROPERTY COMPILE_DEFINITIONS MAKE_PYTHON_MODULE)
set_property(TARGET psi4so PROPERTY LINK_FLAGS "-shared")
set_property(TARGET psi4so PROPERTY OUTPUT_NAME "psi4.so")
target_link_libraries(psi4so ${LINKLIBS})
