find_package(Git QUIET)
if(GIT_FOUND)
   execute_process(COMMAND ${GIT_EXECUTABLE} rev-parse --short HEAD 
                   WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                   OUTPUT_VARIABLE PSI4_HASH
                   OUTPUT_STRIP_TRAILING_WHITESPACE
   )
endif()

set(sources_list version.cc
                 python.cc 
                 psi_start.cc 
                 export_psio.cc 
                 export_mints.cc 
                 psi_stop.cc 
                 export_functional.cc 
                 export_oeprop.cc 
                 export_plugins.cc 
                 export_blas_lapack.cc 
                 export_benchmarks.cc 
                 export_efp.cc 
                 export_cubeprop.cc 
                 clean.cc 
                 create_new_plugin.cc 
                 script.cc 
                 set_memory.cc 
                 read_options.cc 
)
set(PLUGIN_LDFLAGS ${CMAKE_CXX_LINK_FLAGS} ${LIBC_INTERJECT})
add_definitions(-DGIT_VERSION=${PSI4_HASH} 
                -DPSI_VERSION=${PROJECT_VERSION}
                -DPSI_TOP_SRCDIR=${PROJECT_ROOT}
                -Drestrict=${RESTRICT_KEYWORD}
)
set_source_files_properties(create_new_plugin.cc PROPERTIES
                -DPLUGIN_CXX=${CMAKE_CXX_COMPILER}
                -DPLUGIN_INCLUDES=${PSI4_ROOT}
                -DPLUGIN_FLAGS=${CMAKE_CXX_FLAGS}
                -DPLUGIN_DEFINES=""
                -DPLUGIN_OBJDIR=${CMAKE_BINARY_DIR}
                -DPLUGIN_LDFLAGS=${CMAKE_CXX_LINK_FLAGS}
)
psi4_add_binary(psi4_main sources_list mints cubeprop plugin parallel detci
                                       ccresponse ccdensity cchbar fnocc 
                                       cctriples ccenergy cclambda cceom
                                       cctransort ccsort fisapt dfmp2 dcft 
                                       findif mrcc dfocc efp_interface adc
                                       mcscf occ optking
)
