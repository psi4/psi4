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

if(ENABLE_ecpint_RUNTIME)
    if(NOT TARGET ECPINT::ecpint)
        message(FATAL_ERROR "ENABLE_ecpint_RUNTIME requires libecpint package in host environment")
    endif()

    # Include ecpint headers from target WITHOUT creating a dependency (avoid generator expression)
    get_target_property(_ecpint_include_dirs ECPINT::ecpint INTERFACE_INCLUDE_DIRECTORIES)
    target_include_directories(psi4_optional_deps INTERFACE ${_ecpint_include_dirs})
    message(VERBOSE "ECPINT: Using include directories: ${_ecpint_include_dirs}")

    # Link libdl for dlopen/dlsym but do NOT link libecpint itself!
    target_link_libraries(psi4_optional_deps INTERFACE ${CMAKE_DL_LIBS})

    if(UNIX AND NOT APPLE)
        # Calls into libecpint from ecpint.cc (and friends) are left as real
        # *undefined* symbols in core.so; psi4/__init__.py imports core.so
        # with RTLD_LAZY so they're only resolved on first actual ECP use
        # (see libmints/ecpint_loader.cc), not at import time. But many Linux
        # toolchains (notably conda-forge's default, hardened CFLAGS/CXXFLAGS)
        # compile with `-fno-plt`, which makes calls to external functions
        # bypass the PLT entirely and go through a GOT load instead. Those
        # GOT-load relocations are R_X86_64_GLOB_DAT, not R_X86_64_JUMP_SLOT,
        # and GLOB_DAT relocations are *always* resolved eagerly at dlopen()
        # time -- the lazy-binding optimization only applies to PLT/JUMP_SLOT
        # relocations. So under `-fno-plt`, RTLD_LAZY (and even `-Wl,-z,lazy`,
        # see psi4/src/CMakeLists.txt) is silently defeated and `import psi4`
        # crashes with "undefined symbol" again even though libecpint truly
        # isn't installed. GCC/Clang flags are processed left-to-right for
        # toggle options, so appending `-fplt` here (after any inherited
        # `-fno-plt` from CMAKE_CXX_FLAGS) restores real PLT stubs -- and
        # therefore lazy binding -- for calls into libecpint specifically.
        target_compile_options(psi4_optional_deps INTERFACE -fplt)
    endif()

    target_compile_definitions(psi4_optional_deps
      INTERFACE
        USING_ecpint
        USING_ecpint_RUNTIME
      )

elseif(TARGET ECPINT::ecpint)
    target_link_libraries(psi4_optional_deps INTERFACE ECPINT::ecpint)
    target_compile_definitions(psi4_optional_deps
      INTERFACE
        USING_ecpint
        USING_ecpint_STATIC
      )
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
