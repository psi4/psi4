list(APPEND sources
  local.cc
  x2cint.cc
  orbitalspace.cc
  maketab.cc
  oeprop.cc
  writer.cc
  transform.cc
  sieve.cc
  multipolesymmetry.cc
  shellrotation.cc
  deriv.cc
  integralparameters.cc
  twobody.cc
  vector.cc
  sobasis.cc
  cartesianiter.cc
  basisset.cc
  wavefunction.cc
  irrep.cc
  fjt.cc
  chartab.cc
  corrtab.cc
  symop.cc
  benchmark.cc
  get_writer_file_prefix.cc
  petitelist.cc
  solidharmonics.cc
  multipoles.cc
  sointegral.cc
  extern.cc
  factory.cc
  psimath.cc
  dimension.cc
  molecule.cc
  intvector.cc
  mintshelper.cc
  numinthelper.cc
  coordentry.cc
  matrix.cc
  gshell.cc
  integraliter.cc
  pointgrp.cc
  rep.cc
  cdsalclist.cc
  erd_eri.cc
  angularmomentum.cc
  orthog.cc
  )

# l2intf is a listing of all the files that include libint2's boys.h include (which is included in L2's engine.h).  they must
# all be built with the same C++ standard, which is currently chosen as C++14 due to
# problems building with C++17 using Intel compilers on Linux.
# These problems were seen at least through version 2021.4 but were relieved by 2021.8
# Any files including the engine.h should be added here instead of the sources list above.
add_library(l2intf OBJECT 
dipole.cc
electricfield.cc
electrostatic.cc
eribase.cc
eri.cc
integral.cc
kinetic.cc
mcmurchiedavidson.cc  # includes <libint2/boys.h>
multipolepotential.cc # includes <libint2/boys.h>
nabla.cc
onebody.cc
overlap.cc
potential.cc
potentialint.cc
quadrupole.cc
rel_potential.cc
tracelessquadrupole.cc
3coverlap.cc
  )
set_property(TARGET l2intf PROPERTY CXX_STANDARD "${psi4_CXX_STANDARD_Libint2}")
set_property(TARGET l2intf PROPERTY POSITION_INDEPENDENT_CODE ON)
# below are what l2intf files would ordinarily have gotten from psi4_add_module.
# * first brings OpenMP headers
# * second brings l2. Linux needs only headers, but Windows needs headers and linking, so slight overkill here
target_link_libraries(l2intf
  PRIVATE
    tgt::lapack
    Libint2::cxx
  )
target_compile_definitions(l2intf
  PRIVATE
    psi4_SHGSHELL_ORDERING=${psi4_SHGSHELL_ORDERING}
  )

psi4_add_module(lib mints sources)
target_link_libraries(mints PRIVATE l2intf)
target_compile_definitions(mints
  PRIVATE
    psi4_SHGSHELL_ORDERING=${psi4_SHGSHELL_ORDERING}
  )

# Always linked-to external projects
target_link_libraries(mints
  PUBLIC
    Libint2::cxx
  )
# Conditionally linked-to external projects
if(TARGET Libint::libint)
  target_link_libraries(mints
    PUBLIC
      Libint::libint
    )
endif()

if(TARGET ECPINT::ecpint)
  target_sources(mints
    PRIVATE
      ecpint.cc
    )
  target_compile_definitions(mints
    PUBLIC
      USING_ecpint
    )
  target_compile_definitions(l2intf
    PUBLIC
      USING_ecpint
    )
  target_link_libraries(mints
    PUBLIC
      ECPINT::ecpint
    )
  target_link_libraries(l2intf
    PUBLIC
      ECPINT::ecpint
    )
endif()

if(TARGET simint::simint)
  target_sources(l2intf
    PRIVATE
      siminteri.cc
    )
  # Add USING_simint definition, which is not in simint::simint
  target_compile_definitions(l2intf
    PUBLIC
      USING_simint
    )
  target_link_libraries(l2intf
    PUBLIC
      simint::simint
    )
endif()

if(TARGET dkh::dkh)
  target_link_libraries(mints
    PUBLIC
      dkh::dkh
    )
endif()

if(TARGET PCMSolver::pcm)
  target_link_libraries(mints
    PUBLIC
      PCMSolver::pcm
    )
endif()

if(TARGET erd::erd)
  target_link_libraries(mints
    PUBLIC
      erd::erd
    )
endif()

if(TARGET BrianQC::static_wrapper)
  target_compile_definitions(mints
    PUBLIC
      USING_BrianQC
    )
  target_include_directories(mints
    PUBLIC
      $<TARGET_PROPERTY:BrianQC::static_wrapper,INTERFACE_INCLUDE_DIRECTORIES>
    )
endif()
