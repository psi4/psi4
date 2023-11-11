list(APPEND sources
  close_plugin.cc
  load_plugin.cc
  )

psi4_add_module(lib plugin sources)

target_compile_definitions(plugin
  PRIVATE
    $<$<BOOL:${LIBDL_FOUND}>:HAVE_DLFCN_H>
  )

target_link_libraries(plugin
  PUBLIC
    ${LIBDL_LIBRARIES}
  )
