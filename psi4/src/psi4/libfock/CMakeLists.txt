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
target_precompile_headers(fock REUSE_FROM mints)



# Link to common optional dependencies interface library
target_link_libraries(fock PUBLIC psi4_optional_deps)
