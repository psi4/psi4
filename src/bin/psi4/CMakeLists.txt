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

add_executable(psi4 clean.cc create_new_plugin.cc export_benchmarks.cc export_blas_lapack.cc 
                    export_chkpt.cc export_cubefile.cc export_functional.cc export_mints.cc export_oeprop.cc 
                    export_plugins.cc export_psio.cc psi4.cc psi_start.cc psi_stop.cc python.cc
                    read_options.cc script.cc set_memory.cc version.cc)

# Add the python link commands in as raw commands
set(LINKLIBS ${BOOSTLIB} ${PSILIB} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
#set_target_properties(psi4 PROPERTIES LINK_FLAGS ${PYLIB})

#add_dependencies(psi4 ${LINKLIBS})
target_link_libraries(psi4
    ${PSILIB}
    ${Boost_LIBRARIES}
    ${LAPACK_LIBRARIES}
    ${BLAS_LIBRARIES}
    ${PYTHON_LIBRARIES}
)


