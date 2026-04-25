list(APPEND sources
  fisapt.cc
  local2.cc
  )
psi4_add_module(bin fisapt sources)

if(TARGET gauxc::gauxc)
  target_include_directories(fisapt
    PUBLIC
      $<TARGET_PROPERTY:gauxc::gauxc,INTERFACE_INCLUDE_DIRECTORIES>
    )
if (${psi4_ENABLE_PRECOMPILE_HEADERS})
    target_precompile_headers(fisapt REUSE_FROM psio)
endif()

target_link_libraries(fisapt PUBLIC diis)
