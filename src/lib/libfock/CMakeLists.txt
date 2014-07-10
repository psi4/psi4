set(SRC apps.cc cubature.cc direct_screening.cc hamiltonian.cc jk.cc 
        link.cc points.cc solver.cc v.cc wrapper.cc)
if(MPI_FOUND)
   set(SRC ${SRC} Psi4JK.cc)
   message(${PROJECT_SOURCE_DIR}/plugins/JKFactory/src)
   include_directories(${PROJECT_SOURCE_DIR}/plugins/JKFactory/src)
endif()
add_library(fock ${SRC})
add_dependencies(fock mints)
