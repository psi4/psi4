list(APPEND sources
  local.cc
  onebody.cc
  x2cint.cc
  orbitalspace.cc
  osrecur.cc
  maketab.cc
  efpmultipolepotential.cc
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
  eribase.cc
  fjt.cc
  potentialint.cc
  chartab.cc
  corrtab.cc
  quadrupole.cc
  eri.cc
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
  integral.cc
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
  )

psi4_add_module(lib mints sources)

# Always linked-to external projects
target_link_libraries(mints
  PUBLIC
    Libint::libint
  )

# Conditionally linked-to external projects
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
