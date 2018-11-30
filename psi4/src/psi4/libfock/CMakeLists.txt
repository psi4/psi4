list(APPEND sources
  CDJK.cc
  DirectJK.cc
  DiskDFJK.cc
  DiskJK.cc
  GTFockJK.cc
  MemDFJK.cc
  PKJK.cc
  PK_workers.cc
  PKmanagers.cc
  apps.cc
  cubature.cc
  hamiltonian.cc
  jk.cc
  points.cc
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
  )

target_link_libraries(fock
  PRIVATE
    gau2grid::gg
  )
