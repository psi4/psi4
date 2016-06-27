set(sources_list uhf.cc 
                 ks.cc 
                 cuhf.cc 
                 rhf.cc 
                 hf.cc 
                 sad.cc 
                 frac.cc 
                 mom.cc 
                 rohf.cc 
                 stability.cc
)
psi4_add_library(scf_solver sources_list)
target_link_libraries(scf_solver mints)
if(ENABLE_PCMSOLVER)
   target_link_libraries(scf_solver psipcm)
endif()
if(ENABLE_GTFOCK)
  add_dependencies(scf_solver GTFockInt GTFockCInt)
endif()
