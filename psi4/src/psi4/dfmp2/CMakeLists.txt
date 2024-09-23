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
endif()
