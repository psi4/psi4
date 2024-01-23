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
  SplitJK.cc
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

if(TARGET gauxc::gauxc)
  list(APPEND sources
    snLinK.cc
  )
else()
  list(APPEND sources
    snLinK_stub.cc
  )
endif()

psi4_add_module(lib fock sources)

target_compile_definitions(fock
  PRIVATE
    "restrict=${RESTRICT_KEYWORD}"
    $<$<BOOL:${ENABLE_GTFOCK}>:ENABLE_GTFOCK>
    psi4_SHGSHELL_ORDERING=${psi4_SHGSHELL_ORDERING}
  )

if(TARGET gauxc::gauxc)
  target_compile_definitions(fock
    PRIVATE
      "restrict=${RESTRICT_KEYWORD}"
      $<$<BOOL:${ENABLE_GAUXC}>:ENABLE_GAUXC>
    )
endif()

target_link_libraries(fock
  PRIVATE
    gau2grid::gg
    pybind11::headers
    Libint2::cxx  # for <libint2/config.h>
  )

if(TARGET gauxc::gauxc)
  target_include_directories(fock
    PUBLIC
      $<TARGET_PROPERTY:gauxc::gauxc,INTERFACE_INCLUDE_DIRECTORIES>
      #/storage/home/hhive1/dpoole34/data/Codes/GauXC/include/ 
    )
  target_link_libraries(fock
    PUBLIC 
      gauxc::gauxc 
    )
endif()

if(TARGET IntegratorXX::IntegratorXX)
  target_link_libraries(fock
    PRIVATE
      IntegratorXX::IntegratorXX
    )
endif()

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
