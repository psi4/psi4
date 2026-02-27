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
if(ENABLE_GTFOCK)
    add_dependencies(scf_solver INTERFACE GTFockInt GTFockCInt)
endif()

target_link_libraries(scf_solver PUBLIC diis)

if(TARGET ECPINT::ecpint)
  target_compile_definitions(scf_solver
    PUBLIC
      USING_ecpint
    )
  target_link_libraries(scf_solver
    PUBLIC
      ECPINT::ecpint
    )
endif()

if(TARGET gauxc::gauxc)
  target_include_directories(scf_solver
    PUBLIC
      $<TARGET_PROPERTY:gauxc::gauxc,INTERFACE_INCLUDE_DIRECTORIES>
    )
endif()

if(TARGET BrianQC::static_wrapper)
  target_compile_definitions(scf_solver
    PUBLIC
      USING_BrianQC
    )
  target_include_directories(scf_solver
    PUBLIC
      $<TARGET_PROPERTY:BrianQC::static_wrapper,INTERFACE_INCLUDE_DIRECTORIES>
    )
endif()

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
