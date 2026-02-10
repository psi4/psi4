psi4_add_module(bin cc "")

if(TARGET gauxc::gauxc)
  target_include_directories(cc
    PUBLIC
      $<TARGET_PROPERTY:gauxc::gauxc,INTERFACE_INCLUDE_DIRECTORIES>
    )
endif()

# FIXME Remove once done
list(APPEND _Clang_debug_flags
  "-Wall"
  "-O0"
  "-g3"
  "-Wextra"
  "-Winit-self"
  "-Woverloaded-virtual"
  "-Wunused"
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
target_sources(cc
  PRIVATE
  amps.cc
  )
target_compile_options(cc
  PRIVATE
  $<$<AND:$<CXX_COMPILER_ID:Clang>,$<CONFIG:Debug>>:${_Clang_debug_flags}>
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
