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
  matrix_eigen.cc
  gshell.cc
  integraliter.cc
  pointgrp.cc
  rep.cc
  cdsalclist.cc
  angularmomentum.cc
  orthog.cc
  thc_eri.cc
  )

# l2intf is a listing of all the files that include libint2's boys.h include (which is included in L2's engine.h).
# These were originally separated so that a separate C++ standard could be applied, as Intel compilers were having problems at least through version 2021.4 (relieved by 2021.8)
# If that ever needs to be revived, look for psi4_CXX_STANDARD_Libint2 in history.
# Nowadays, these are files that take a lot of memory to compile and are prone to silently failing on CI,
#   esp. with Ninja + Linux + CI-or-laptop, so they're segregated to throttle build concurrence.
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
set_property(TARGET l2intf PROPERTY POSITION_INDEPENDENT_CODE ON)

# below are what l2intf files would ordinarily have gotten from psi4_add_module.
target_link_libraries(l2intf
  PRIVATE
    psi4_required_deps
    psi4_optional_deps
  )

psi4_add_module(lib mints sources)
target_link_libraries(mints PRIVATE l2intf)
add_dependencies(mints l2intf)  # impose build-order dependence allowing CMAKE_OPTIMIZE_DEPENDENCIES=ON

if (psi4_ENABLE_PRECOMPILE_HEADERS)
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

# Limit l2intf compile memory to avoid OOM, esp. Ninja on CI/cloud runners.
# Measured (GCC 14, -O3, L2 2.14): every l2intf file that includes libint2/engine.h peaks at 5.7-7.5 GB
#   by itself, while one unity TU of all 19 files peaks at 7.6 GB in ~200 s -- no more memory or time than
#   a single heavy file, vs. 3 x 7.5 GB (and ~3x the CPU) for 8-file batches.
#   clang 20 has the same pattern at ~4.0-4.4 GB per TU (single file or all 19), with similar times.
# (1) Compile l2intf as unity TUs of psi4_l2intf_UNITY_BUILD_BATCH_SIZE files, independent of CMAKE_UNITY_BUILD
#   and CMAKE_UNITY_BUILD_BATCH_SIZE. Default 0 = no limit, i.e. one TU for all files: least memory and CPU,
#   right for CI and one-shot builds, but editing any l2intf file recompiles all of them (~3.5 min).
#   Developers editing these files can set 1 (per-file TUs; ~1 min to rebuild one file, mostly engine.h).
if(NOT DEFINED psi4_l2intf_UNITY_BUILD_BATCH_SIZE)
    set(psi4_l2intf_UNITY_BUILD_BATCH_SIZE 0)
endif()
if(NOT psi4_l2intf_UNITY_BUILD_BATCH_SIZE MATCHES "^[0-9]+$")
    message(FATAL_ERROR "psi4_l2intf_UNITY_BUILD_BATCH_SIZE must be a non-negative integer, not '${psi4_l2intf_UNITY_BUILD_BATCH_SIZE}'")
endif()
set_property(TARGET l2intf PROPERTY UNITY_BUILD ON)
set_property(TARGET l2intf PROPERTY UNITY_BUILD_BATCH_SIZE ${psi4_l2intf_UNITY_BUILD_BATCH_SIZE})
message(STATUS "l2intf unity batch size: ${psi4_l2intf_UNITY_BUILD_BATCH_SIZE} (0 = one TU)")

# (2) With unity build off (or Makefile generators, which ignore pools), each l2intf file is its own 6-8 GB job.
#   For Ninja, cap concurrent l2intf compiles by memory: one job per ~8 GB available, at least one.
#   Use the container's cgroup limit if it's tighter than the host memory CMake reports.
#   The pool is shared by all compiles that include libint2/engine.h: l2intf here, fock_l2engine
#   (libfock/snLinK.cc), and the pybind11 core module (core.cc, export_mints.cc; psi4/src/CMakeLists.txt).
cmake_host_system_information(RESULT _l2intf_mem QUERY AVAILABLE_PHYSICAL_MEMORY)  # MiB
foreach(_cg_max _cg_cur IN ZIP_LISTS
        "/sys/fs/cgroup/memory.max;/sys/fs/cgroup/memory/memory.limit_in_bytes"
        "/sys/fs/cgroup/memory.current;/sys/fs/cgroup/memory/memory.usage_in_bytes")
    if(EXISTS "${_cg_max}" AND EXISTS "${_cg_cur}")
        file(READ "${_cg_max}" _cg_max_val)
        file(READ "${_cg_cur}" _cg_cur_val)
        string(STRIP "${_cg_max_val}" _cg_max_val)
        string(STRIP "${_cg_cur_val}" _cg_cur_val)
        # "max" (v2) means unlimited; v1 unlimited is a huge number that the comparison below discards
        if(_cg_max_val MATCHES "^[0-9]+$" AND _cg_cur_val MATCHES "^[0-9]+$")
            math(EXPR _cg_avail "(${_cg_max_val} - ${_cg_cur_val}) / 1048576")
            if(_cg_avail LESS _l2intf_mem)
                set(_l2intf_mem ${_cg_avail})
            endif()
        endif()
        break()
    endif()
endforeach()
math(EXPR _l2intf_jobs "${_l2intf_mem} / 8192")
if(_l2intf_jobs LESS 1)
    set(_l2intf_jobs 1)
endif()
set_property(GLOBAL APPEND PROPERTY JOB_POOLS libint2_engine_pool=${_l2intf_jobs})
set_property(TARGET l2intf PROPERTY JOB_POOL_COMPILE libint2_engine_pool)
message(STATUS "libint2_engine_pool: ${_l2intf_mem} MiB available -> compile job pool depth ${_l2intf_jobs}")

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

if(TARGET cint::cint)
  target_sources(mints
    PRIVATE
      libcinteri.cc
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
