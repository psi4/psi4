set(sources_list apps.cc
                 v.cc
                 hamiltonian.cc
                 points.cc
                 cubature.cc
                 solver.cc
                 wrapper.cc
                 jk.cc
                 DiskJK.cc
                 PKJK.cc
                 DirectJK.cc
                 DFJK.cc
                 CDJK.cc
                 GTFockJK.cc
                 soscf.cc
                 PK_workers.cc
                 PKmanagers.cc
                 symm_JK.cc
)
add_definitions("-Drestrict=${RESTRICT_KEYWORD}")

if(ENABLE_GTFOCK)
   add_definitions("-DENABLE_GTFOCK")
endif()

psi4_add_module(lib fock sources_list mints functional 3index psio)
target_link_libraries(fock PRIVATE gau2grid::gg)
