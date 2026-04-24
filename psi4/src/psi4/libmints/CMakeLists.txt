list(APPEND sources
  local.cc
  x2cint.cc
  zora.cc
  orbitalspace.cc
  maketab.cc
  oeprop.cc
  writer.cc
  transform.cc
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
  angularmomentum.cc
  orthog.cc
  thc_eri.cc
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
potential_erf.cc
quadrupole.cc
rel_potential.cc
tracelessquadrupole.cc
3coverlap.cc
  )
#set_property(TARGET l2intf PROPERTY CXX_STANDARD "${psi4_CXX_STANDARD_Libint2}")
set_property(TARGET l2intf PROPERTY POSITION_INDEPENDENT_CODE ON)

# below are what l2intf files would ordinarily have gotten from psi4_add_module.
# l2intf is an OBJECT library, so it doesn't automatically get the dependencies from the macro.
# We need to explicitly link the required and optional deps here.
# Note: OBJECT libraries need PRIVATE for their own compilation but will still propagate
# usage requirements to targets that link against them.
target_link_libraries(l2intf
  PRIVATE
    psi4_required_deps
    psi4_optional_deps
  )

psi4_add_module(lib mints sources)
target_link_libraries(mints PRIVATE l2intf)
add_dependencies(mints l2intf)  # allow CMAKE_OPTIMIZE_DEPENDENCIES=ON
target_compile_definitions(mints
  PRIVATE
    psi4_SHGSHELL_ORDERING=${psi4_SHGSHELL_ORDERING}
  )

if (${psi4_ENABLE_PRECOMPILE_HEADERS})
   # When OpenOrbitalOptimizer is enabled, mints cannot reuse psio's PCH because:
   # 1. matrix.h conditionally includes <armadillo> based on USING_OpenOrbitalOptimizer
   # 2. psio's PCH sees matrix.h without USING_OpenOrbitalOptimizer (only forward decls)
   # 3. mints sees matrix.h WITH USING_OpenOrbitalOptimizer (full armadillo)
   # 4. This creates typedef conflicts and incomplete type errors
   # Solution: Skip PCH for mints when OpenOrbitalOptimizer is enabled
   # Note: Adding USING_OpenOrbitalOptimizer to psio doesn't work because it causes
   # BLAS symbol conflicts in unrelated modules (fnocc, libqt) that also use psio's PCH
   if(NOT TARGET OpenOrbitalOptimizer::OpenOrbitalOptimizer)
      target_precompile_headers(mints REUSE_FROM psio)
      target_precompile_headers(l2intf REUSE_FROM psio)
   endif()
endif()

# limit to avoid OOM errors with Ninja parallelism
# * finer control than `CMAKE_BUILD_PARALLEL_LEVEL=1 cmake --build`
cmake_host_system_information(RESULT _memfree QUERY AVAILABLE_PHYSICAL_MEMORY)
cmake_host_system_information(RESULT _numcore QUERY NUMBER_OF_PHYSICAL_CORES)
math(EXPR _memrat "${_memfree} / ${_numcore}")
message("AVAILABLE_PHYSICAL_MEMORY ${_memfree} / NUMBER_OF_PHYSICAL_CORES ${_numcore} = ${_memrat}")
if(_numcore LESS 5)
    if (UNIX AND NOT APPLE)
        set_property(TARGET l2intf PROPERTY JOB_POOL_COMPILE one_job)
    endif()
else()
    set_property(TARGET l2intf PROPERTY JOB_POOL_COMPILE four_jobs)
endif()

# Conditionally add sources that depend on optional dependencies
if(TARGET ECPINT::ecpint)
  target_sources(mints
    PRIVATE
      ecpint.cc
    )
endif()

if(TARGET simint::simint)
  target_sources(l2intf
    PRIVATE
      siminteri.cc
    )
endif()

# OpenOrbitalOptimizer needs to be PUBLIC so USING_OpenOrbitalOptimizer propagates
# to modules that link to libmints (this enables armadillo support in matrix.h)
# but doesn't pollute unrelated modules like libqt (which would get BLAS conflicts)
if(TARGET OpenOrbitalOptimizer::OpenOrbitalOptimizer)
  target_compile_definitions(mints
    PUBLIC
      USING_OpenOrbitalOptimizer
    )
endif()
