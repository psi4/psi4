list(APPEND sources
  brianqc_int.cc
  CDJK.cc
  CompositeJK.cc
  COSK.cc
  DirectDFJ.cc
  DirectJK.cc
  DiskDFJK.cc
  DiskJK.cc
  gauxc_int.cc
  GTFockJK.cc
  integrator_manager.cc
  LinK.cc
  MemDFJK.cc
  PKJK.cc
  PK_workers.cc
  PKmanagers.cc
  SplitJK.cc
  apps.cc
  cubature.cc
  hamiltonian.cc
  jk.cc
  points.cc
  sap.cc
  snLinK.cc
  solver.cc
  soscf.cc
  v.cc
  wrapper.cc
  )

psi4_add_module(lib fock sources)

target_compile_definitions(fock
  PRIVATE
    "restrict=${RESTRICT_KEYWORD}"
    $<$<BOOL:${ENABLE_GTFOCK}>:ENABLE_GTFOCK>
    psi4_SHGSHELL_ORDERING=${psi4_SHGSHELL_ORDERING}
  )

target_link_libraries(fock
  PRIVATE
    gau2grid::gg
    pybind11::headers
    Libint2::cxx  # for <libint2/config.h>
  )

# Link to common optional dependencies interface library
target_link_libraries(fock PUBLIC psi4_optional_deps)
