list(APPEND sources
  cuhf.cc
  frac.cc
  hf.cc
  mom.cc
  rhf.cc
  rohf.cc
  sad.cc
  uhf.cc
  cghf.cc
  )


psi4_add_module(lib scf_solver sources)
if (psi4_ENABLE_PRECOMPILE_HEADERS)
    # Skip PCH when OpenOrbitalOptimizer is enabled to avoid conflicts with armadillo
    if(NOT TARGET OpenOrbitalOptimizer::OpenOrbitalOptimizer)
        target_precompile_headers(scf_solver REUSE_FROM psio)
    endif()
endif()
if(ENABLE_GTFOCK)
    add_dependencies(scf_solver INTERFACE GTFockInt GTFockCInt)
endif()

target_link_libraries(scf_solver PUBLIC diis)

# OpenOrbitalOptimizer support for SCF solver
# PUBLIC because scf_solver headers may expose OpenOrbitalOptimizer types
if(TARGET OpenOrbitalOptimizer::OpenOrbitalOptimizer)
  target_compile_definitions(scf_solver
    PUBLIC
      USING_OpenOrbitalOptimizer
    )
  target_link_libraries(scf_solver
    PUBLIC
      OpenOrbitalOptimizer::OpenOrbitalOptimizer
    )
endif()

if(TARGET Einsums::Einsums)
  target_compile_definitions(scf_solver
    PUBLIC
    USING_Einsums
  )
  target_link_libraries(scf_solver PUBLIC Einsums::Einsums)

  if(DEFINED ENV{CONDA_TOOLCHAIN_BUILD})
    if(APPLE)
      # * possibly these should be exported from Einsums somehow
      # * when using conda clang, compiling against an old SDK &
      #   an openmp symbols is absent
      # * afaik, envvar above is defined when `conda activate envvar_w_compilers`.
      #   There's a CONDA_BACKUP_CONDA_TOOLCHAIN_BUILD envvar when not activated.
      #   So hopefully, this doesn't interfere with system compilers.
      target_compile_definitions(
        scf_solver
        PUBLIC
          _LIBCPP_DISABLE_AVAILABILITY
          H5CPP_USE_OMP_ALIGNED_ALLOC
        )
  endif()
endif()

endif()
