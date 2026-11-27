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

psi4_add_module(lib fock sources)
if (psi4_ENABLE_PRECOMPILE_HEADERS)
    target_precompile_headers(fock REUSE_FROM psio)
endif()

# example for turning off unity build for a directory while developint
# set_property(TARGET fock PROPERTY UNITY_BUILD OFF)

target_link_libraries(fock
  PRIVATE
    pybind11::headers
  )

# snLinK.cc includes <libint2.hpp>: with libint2 >=2.14 such a TU peaks at 5.6+ GB (GCC) just from the
#   libint2 engine headers. Compile it separately so it shares the memory-sized job pool with the
#   libmints l2intf files (see libmints/CMakeLists.txt) instead of joining a fock unity batch.
add_library(fock_l2engine OBJECT snLinK.cc)
set_target_properties(fock_l2engine
  PROPERTIES
    POSITION_INDEPENDENT_CODE ON
    CXX_VISIBILITY_PRESET hidden
    VISIBILITY_INLINES_HIDDEN 1
    JOB_POOL_COMPILE libint2_engine_pool
  )
target_link_libraries(fock_l2engine
  PRIVATE
    psi4_required_deps
    psi4_optional_deps
    pybind11::headers
  )
if (psi4_ENABLE_PRECOMPILE_HEADERS)
    target_precompile_headers(fock_l2engine REUSE_FROM psio)
endif()
target_link_libraries(fock PRIVATE fock_l2engine)
add_dependencies(fock fock_l2engine)  # impose build-order dependence allowing CMAKE_OPTIMIZE_DEPENDENCIES=ON (Makefiles need it)
