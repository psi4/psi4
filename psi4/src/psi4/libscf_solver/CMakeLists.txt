set(sources_list uhf.cc 
                 cuhf.cc 
                 rhf.cc 
                 hf.cc 
                 sad.cc 
                 frac.cc 
                 mom.cc 
                 rohf.cc 
                 stability.cc
)
psi4_add_module(lib scf_solver sources_list mints diis efp_solver fock)
if(TARGET PCMSolver::PCMSolver)
    # external libpcm needs host_writer from psipcm, hence the forced symbol export
    target_link_libraries(scf_solver PRIVATE ${PRE_LIBRARY_OPTION} psipcm ${POST_LIBRARY_OPTION})
    target_compile_definitions(scf_solver PRIVATE $<TARGET_PROPERTY:PCMSolver::PCMSolver,INTERFACE_COMPILE_DEFINITIONS>)
endif()
if(ENABLE_GTFOCK)
    add_dependencies(scf_solver INTERFACE GTFockInt GTFockCInt)
endif()
