add_executable(psi4 clean.cc create_new_plugin.cc export_benchmarks.cc export_blas_lapack.cc 
                    export_chkpt.cc export_functional.cc export_mints.cc export_oeprop.cc 
                    export_plugins.cc export_psio.cc psi4.cc psi_start.cc psi_stop.cc python.cc
                    read_options.cc script.cc set_memory.cc version.cc)
message(${BOOSTLIB})
target_link_libraries(psi4 ${BOOSTLIB} ${PYLIB})
include_directories(${BOOSTINC} ${PYINC})
