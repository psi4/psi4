list(APPEND sources
  mp2.cc
  corr_grad.cc
  wrapper.cc
  )
psi4_add_module(bin dfmp2 sources)

if(TARGET gauxc::gauxc)
  target_include_directories(dfmp2
    PUBLIC
      $<TARGET_PROPERTY:gauxc::gauxc,INTERFACE_INCLUDE_DIRECTORIES>
    )
if (${psi4_ENABLE_PRECOMPILE_HEADERS})
    target_precompile_headers(dfmp2 REUSE_FROM psio)
endif()
