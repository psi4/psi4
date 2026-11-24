set(sources
  dftensor.cc
  dfhelper.cc
  denominator.cc
  fittingmetric.cc
  cholesky.cc
  )

if(TARGET cuEST::cuEST)
    list(APPEND sources
        cuESTMOInts.cc
        )
endif()

psi4_add_module(lib 3index sources)

if (psi4_ENABLE_PRECOMPILE_HEADERS)
    target_precompile_headers(3index REUSE_FROM psio)
endif()
