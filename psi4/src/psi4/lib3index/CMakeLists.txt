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
if (${psi4_ENABLE_PRECOMPILE_HEADERS})
    target_precompile_headers(3index REUSE_FROM psio)
endif()
