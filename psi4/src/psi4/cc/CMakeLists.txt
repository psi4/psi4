add_library(cc STATIC
  ccwave.cc
  ccmoinfo.cc
  ccparams.cc
  )

# FIXME Remove once done
list(APPEND _Clang_debug_flags
  "-Wall"
  "-O0"
  "-g3"
  "-Wextra"
  "-Winit-self"
  "-Woverloaded-virtual"
  "-Wuninitialized"
  "-Wmissing-declarations"
  "-Wwrite-strings"
  "-Weffc++"
  "-Wdocumentation"
  "-Wno-sign-compare"
  "-m64"
  "-Qunused-arguments"
  "-fcolor-diagnostics"
  )
target_compile_options(cc
  PRIVATE
  $<$<AND:$<CXX_COMPILER_ID:Clang>,$<CONFIG:Debug>>:${_Clang_debug_flags}>
  )

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
