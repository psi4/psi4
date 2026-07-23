# Create interface library for required dependencies
add_library(psi4_required_deps INTERFACE)

target_link_libraries(psi4_required_deps
  INTERFACE
    tgt::lapack     # Provides OpenMP flags and BLAS/LAPACK headers
    Libint2::cxx    # Required for integral engine
    Libxc::xc       # Required for DFT functionals; target provides USING_Libxc
    gau2grid::gg    # Required for DFT grid generation
)

set (CMAKE_OPTIMIZE_DEPENDENCIES ON)

# Create interface library for common optional dependencies to avoid repeating if(TARGET...) blocks
add_library(psi4_optional_deps INTERFACE)

if(TARGET gauxc::gauxc)
  target_link_libraries(psi4_optional_deps INTERFACE gauxc::gauxc)
  target_compile_definitions(psi4_optional_deps INTERFACE USING_gauxc)
endif()

if(TARGET ECPINT::ecpint)
  target_link_libraries(psi4_optional_deps INTERFACE ECPINT::ecpint)
  target_compile_definitions(psi4_optional_deps INTERFACE USING_ecpint)
endif()

if(TARGET OpenOrbitalOptimizer::OpenOrbitalOptimizer)
  # Note: USING_OpenOrbitalOptimizer is NOT added to psi4_optional_deps because:
  # 1. USING_OpenOrbitalOptimizer causes armadillo headers to be included in matrix.h
  # 2. Armadillo declares BLAS functions which conflict with libqt's BLAS declarations
  # 3. Only libmints, libscf_solver, scfgrad actually need OpenOrbitalOptimizer
  # Solution: libmints adds USING_ as PUBLIC, which keeps it away from upstream modules like libqt
  target_link_libraries(psi4_optional_deps INTERFACE OpenOrbitalOptimizer::OpenOrbitalOptimizer)
endif()

if(TARGET BrianQC::static_wrapper)
  target_link_libraries(psi4_optional_deps INTERFACE BrianQC::static_wrapper)
  target_compile_definitions(psi4_optional_deps INTERFACE USING_BrianQC)
endif()

if(TARGET IntegratorXX::IntegratorXX)
  target_link_libraries(psi4_optional_deps INTERFACE IntegratorXX::IntegratorXX)
  target_compile_definitions(psi4_optional_deps INTERFACE USING_IntegratorXX)
endif()

if(TARGET simint::simint)
  target_link_libraries(psi4_optional_deps INTERFACE simint::simint)
  target_compile_definitions(psi4_optional_deps INTERFACE USING_simint)
endif()

if(TARGET libcint::libcint)
  target_link_libraries(psi4_optional_deps INTERFACE libcint::libcint)
  target_compile_definitions(psi4_optional_deps INTERFACE USING_libcint)
endif()

if(TARGET dkh::dkh)
  target_link_libraries(psi4_optional_deps INTERFACE dkh::dkh)
  # USING_dkh provided by dkh::dkh target itself
endif()

if(TARGET PCMSolver::pcm)
  target_link_libraries(psi4_optional_deps INTERFACE PCMSolver::pcm)
  # USING_PCMSolver provided by PCMSolver::pcm target itself
endif()

if(TARGET CheMPS2::chemps2)
  target_link_libraries(psi4_optional_deps INTERFACE CheMPS2::chemps2)
  # USING_CheMPS2 provided by CheMPS2::chemps2 target itself
endif()

if(TARGET ambit::ambit)
  target_link_libraries(psi4_optional_deps INTERFACE ambit::ambit)
  # USING_ambit provided by ambit::ambit target itself
endif()

if(TARGET Einsums::Einsums)
  target_link_libraries(psi4_optional_deps INTERFACE Einsums::Einsums)
  target_compile_definitions(psi4_optional_deps INTERFACE USING_Einsums)
endif()

foreach(dir_name
  lib3index
  libciomr
  libcubeprop
  libdiis
  libdisp
  libdpd
  libfilesystem
  libfock
  libfunctional
  libiwl
  libmints
  libmoinfo
  liboptions
  libplugin
  libpsi4util
  libpsio
  libpsipcm
  libqt
  libsapt_solver
  libscf_solver
  libtrans
  )
  add_subdirectory(${dir_name})
endforeach ()

foreach(dir_name
  ambit_interface
  cc
  cctransort
  dct
  detci
  dfep2
  dfmp2
  dfocc
  dlpno
  dmrg
  fisapt
  fnocc
  mcscf
  mrcc
  occ
  psimrcc
  sapt
  scfgrad
  f12
  dummy_einsums
  dummy_integratorxx
  )
  add_subdirectory(${dir_name})
endforeach ()
