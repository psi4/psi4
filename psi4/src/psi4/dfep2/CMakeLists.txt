list(APPEND sources
  dfep2.cc
  )
psi4_add_module(bin dfep2 sources)

if(TARGET gauxc::gauxc)
  target_include_directories(dfep2
    PUBLIC
      $<TARGET_PROPERTY:gauxc::gauxc,INTERFACE_INCLUDE_DIRECTORIES>
    )
endif()
