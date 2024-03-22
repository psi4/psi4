set(sources
  dftensor.cc
  dfhelper.cc
  denominator.cc
  fittingmetric.cc
  cholesky.cc
  )
psi4_add_module(lib 3index sources)

if(TARGET gauxc::gauxc)
  target_include_directories(3index
    PUBLIC
      $<TARGET_PROPERTY:gauxc::gauxc,INTERFACE_INCLUDE_DIRECTORIES>
    )
  target_link_libraries(3index
    PUBLIC 
      gauxc::gauxc 
    )
endif()
