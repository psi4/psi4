#
# Run the git version script, and update gitversion.h, if it's changed
#
add_custom_command(
    COMMAND ${PYTHON_EXECUTABLE} gitversion.py
    COMMAND ${CMAKE} -E copy_if_different gitversion.h.tmp gitversion.h
    OUTPUT gitversion.h
    COMMENT "Generating Git info"
)

add_executable(psi4 clean.cc create_new_plugin.cc export_benchmarks.cc export_blas_lapack.cc 
                    export_chkpt.cc export_functional.cc export_mints.cc export_oeprop.cc 
                    export_plugins.cc export_psio.cc psi4.cc psi_start.cc psi_stop.cc python.cc
                    read_options.cc script.cc set_memory.cc version.cc)

target_link_libraries(psi4 ${BOOSTLIB} ${PYLIB} ${PSILIB})
