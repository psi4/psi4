list(APPEND sources
  jk_grad.cc
  response.cc
  scf_grad.cc
  wrapper.cc
  )
psi4_add_module(bin scfgrad sources)

if(TARGET BrianQC::static_wrapper)
  target_compile_definitions(scfgrad
    PUBLIC
      USING_BrianQC
    )
  target_include_directories(scfgrad
    PUBLIC
      $<TARGET_PROPERTY:BrianQC::static_wrapper,INTERFACE_INCLUDE_DIRECTORIES>
    )
endif()
