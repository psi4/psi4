list(APPEND sources
  cuhf.cc
  frac.cc
  hf.cc
  mom.cc
  rhf.cc
  rohf.cc
  sad.cc
  uhf.cc
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
