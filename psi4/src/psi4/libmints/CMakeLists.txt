list(APPEND sources
  local.cc
  onebody.cc
  x2cint.cc
  orbitalspace.cc
  osrecur.cc
  maketab.cc
  multipolepotential.cc
  rel_potential.cc
  oeprop.cc
  writer.cc
  transform.cc
  sieve.cc
  multipolesymmetry.cc
  shellrotation.cc
  deriv.cc
  overlap.cc
  integralparameters.cc
  twobody.cc
  vector.cc
  sobasis.cc
  cartesianiter.cc
  basisset.cc
  electrostatic.cc
  wavefunction.cc
  irrep.cc
  fjt.cc
  potentialint.cc
  chartab.cc
  corrtab.cc
  quadrupole.cc
  symop.cc
  benchmark.cc
  get_writer_file_prefix.cc
  3coverlap.cc
  petitelist.cc
  solidharmonics.cc
  electricfield.cc
  multipoles.cc
  dipole.cc
  sointegral.cc
  extern.cc
  nabla.cc
  factory.cc
  psimath.cc
  dimension.cc
  molecule.cc
  intvector.cc
  potential.cc
  mintshelper.cc
  coordentry.cc
  kinetic.cc
  tracelessquadrupole.cc
  pseudospectral.cc
  matrix.cc
  gshell.cc
  integraliter.cc
  pointgrp.cc
  rep.cc
  cdsalclist.cc
  erd_eri.cc
  angularmomentum.cc
  bessel.cc
  gaussquad.cc
  ecpint.cc
  orthog.cc
  eri.cc      # ~TMP
  eribase.cc  # ~TMP
  integral.cc # ~TMP
  )

if (MSVC)
    set_source_files_properties(eri.cc eribase.cc integral.cc PROPERTIES COMPILE_OPTIONS "-std:c++14")
else()
    set_source_files_properties(eri.cc eribase.cc integral.cc PROPERTIES COMPILE_OPTIONS "-std=c++14")
endif()

#TMP add_library(l2intf OBJECT eri.cc eribase.cc integral.cc)
#TMP set_property(TARGET l2intf PROPERTY CXX_STANDARD 14)
#TMP set_property(TARGET l2intf PROPERTY POSITION_INDEPENDENT_CODE ON)
#TMP target_include_directories(l2intf
#TMP   PRIVATE
#TMP     $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
#TMP     $<TARGET_PROPERTY:Libint2::cxx,INTERFACE_INCLUDE_DIRECTORIES>
#TMP   )

psi4_add_module(lib mints sources)
#TMP target_link_libraries(mints PRIVATE l2intf)

# Always linked-to external projects
#target_link_libraries(l2intf
#  PUBLIC
#    Libint2::cxx
#  )
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

if(TARGET simint::simint)
  target_sources(mints
    PRIVATE
      siminteri.cc
    )
  # Add USING_simint definition, which is not in simint::simint
  target_compile_definitions(mints
    PUBLIC
      USING_simint
    )
  target_link_libraries(mints
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

if(TARGET ECPINT::ecpint)
  target_link_libraries(mints
    PUBLIC
      ECPINT::ecpint
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
