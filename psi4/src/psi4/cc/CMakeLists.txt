add_library(cc STATIC "")

set_target_properties(cc
  PROPERTIES
    POSITION_INDEPENDENT_CODE 1
  )

list(APPEND _ccmodules
  ccdensity
  ccenergy
  cceom
  cchbar
  cclambda
  ccresponse
  cctriples
  )
foreach(dir_name IN LISTS _ccmodules)
  add_subdirectory(${dir_name})
endforeach ()

set_property(GLOBAL APPEND PROPERTY BINLIST cc)

target_link_libraries(cc
  PRIVATE
    mints
    tgt::lapack
  )
