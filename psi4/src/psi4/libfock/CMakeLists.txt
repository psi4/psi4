list(APPEND sources
  CDJK.cc
  CompositeJK.cc
  COSK.cc
  DirectDFJ.cc
  DirectJK.cc
  DiskDFJK.cc
  DiskJK.cc
  GTFockJK.cc
  LinK.cc
  MemDFJK.cc
  PKJK.cc
  PK_workers.cc
  PKmanagers.cc
  apps.cc
  cubature.cc
  hamiltonian.cc
  jk.cc
  points.cc
  sap.cc
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

if(TARGET BrianQC::static_wrapper)
  target_compile_definitions(fock
    PUBLIC
      USING_BrianQC
    )
  target_include_directories(fock
    PUBLIC
      $<TARGET_PROPERTY:BrianQC::static_wrapper,INTERFACE_INCLUDE_DIRECTORIES>
    )
endif()
