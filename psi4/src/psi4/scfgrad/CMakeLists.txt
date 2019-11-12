list(APPEND sources
  jk_grad.cc
  response.cc
  scf_grad.cc
  wrapper.cc
  )
psi4_add_module(bin scfgrad sources)

if(ENABLE_BrianQC)
  target_compile_definitions(scfgrad
    PUBLIC
      USING_BrianQC
    )
endif()
