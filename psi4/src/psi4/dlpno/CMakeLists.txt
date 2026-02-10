list(APPEND sources
  dlpno.cc
  mp2.cc
  ccsd.cc
  triples.cc
  wrapper.cc
  sparse.cc
  )
psi4_add_module(bin dlpno sources)

if(TARGET gauxc::gauxc)
  target_include_directories(dlpno
    PUBLIC
      $<TARGET_PROPERTY:gauxc::gauxc,INTERFACE_INCLUDE_DIRECTORIES>
    )
endif()

target_link_libraries(dlpno PUBLIC diis)
