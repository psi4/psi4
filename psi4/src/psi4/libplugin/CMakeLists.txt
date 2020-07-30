list(APPEND sources
  close_plugin.cc
  load_plugin.cc
  )

psi4_add_module(lib plugin sources)

if(CMAKE_DL_LIBS)
  target_compile_definitions(plugin
    PRIVATE
      HAVE_DLFCN_H
    )

  target_link_libraries(plugin
    PUBLIC
      ${CMAKE_DL_LIBS}
    )
endif()
