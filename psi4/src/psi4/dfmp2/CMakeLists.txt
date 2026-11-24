list(APPEND sources
  mp2.cc
  corr_grad.cc
  wrapper.cc
  )

if(TARGET cuEST::cuEST)
    list(APPEND sources 
        mp2_kernels.cu
        )
endif()

psi4_add_module(bin dfmp2 sources)

if (psi4_ENABLE_PRECOMPILE_HEADERS)
    target_precompile_headers(dfmp2 REUSE_FROM psio)
endif()

