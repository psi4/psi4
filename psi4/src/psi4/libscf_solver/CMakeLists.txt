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
if(TARGET PCMSolver::pcm)
    target_link_libraries(scf_solver PRIVATE psipcm PCMSolver::pcm)
endif()
if(ENABLE_GTFOCK)
    add_dependencies(scf_solver INTERFACE GTFockInt GTFockCInt)
endif()
