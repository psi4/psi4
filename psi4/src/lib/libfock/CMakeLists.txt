set(sources_list apps.cc 
                 v.cc 
                 hamiltonian.cc 
                 points.cc
                 cubature.cc 
                 solver.cc 
                 link.cc 
                 direct_screening.cc
                 wrapper.cc 
                 jk.cc 
                 DiskJK.cc 
                 PKJK.cc 
                 DirectJK.cc 
                 DFJK.cc
                 CDJK.cc 
                 FastDFJK.cc 
                 PSJK.cc 
                 GTFockJK.cc 
                 soscf.cc
                 PK_workers.cc
                 PKmanagers.cc
)
add_definitions("-Drestrict=${RESTRICT_KEYWORD}")
psi4_add_library(fock sources_list mints functional 3index psio)
target_include_directories(fock PUBLIC ${Boost_INCLUDE_DIRS})
