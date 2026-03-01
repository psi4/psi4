list(APPEND sources
  csg.cc
  cubeprop.cc
  )
psi4_add_module(lib cubeprop sources)

if(TARGET cuEST::cuEST)
  target_link_libraries(cubeprop
    PUBLIC
      cuEST::cuEST
    )
endif()

