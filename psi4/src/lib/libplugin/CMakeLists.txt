set(sources_list close_plugin.cc 
                 load_plugin.cc
)
if(LIBDL_FOUND)
   add_definitions(-DHAVE_DLFCN_H)
endif()
psi4_add_library(plugin sources_list mints ${LIBDL_LIBRARIES})