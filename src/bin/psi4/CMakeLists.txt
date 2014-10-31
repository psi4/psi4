#
# Run the git version script, and update gitversion.h, if it's changed
#
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

set(SRC clean.cc create_new_plugin.cc export_benchmarks.cc export_blas_lapack.cc export_libfrag.cc export_chkpt.cc export_cubefile.cc export_functional.cc export_mints.cc export_oeprop.cc export_plugins.cc export_psio.cc psi4.cc psi_start.cc psi_stop.cc python.cc read_options.cc script.cc set_memory.cc version.cc)

set(LINKLIBS ${PSILIB} ${Boost_LIBRARIES} ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} ${PYTHON_LIBRARIES} dl)

# executable psi4
add_executable(psi4 ${SRC})
if(CUSTOM_BOOST_BUILD)
   add_dependencies(psi4 custom_boost)
endif()
target_link_libraries(psi4 ${LINKLIBS})

# standalone python module psi4.so
add_executable(psi4so EXCLUDE_FROM_ALL ${SRC})
set_property(TARGET psi4so PROPERTY COMPILE_DEFINITIONS MAKE_PYTHON_MODULE)
set_property(TARGET psi4so PROPERTY LINK_FLAGS "-shared")
set_property(TARGET psi4so PROPERTY OUTPUT_NAME "psi4.so")
target_link_libraries(psi4so ${LINKLIBS})
