list(APPEND sources
  csg.cc
  cubeprop.cc
  )
psi4_add_module(lib cubeprop sources)

if(TARGET gauxc::gauxc)
  target_include_directories(cubeprop
    PUBLIC
      $<TARGET_PROPERTY:gauxc::gauxc,INTERFACE_INCLUDE_DIRECTORIES>
    )
endif()
