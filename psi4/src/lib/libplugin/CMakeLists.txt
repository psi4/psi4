set(sources_list close_plugin.cc 
                 load_plugin.cc
)
psi4_add_library(plugin sources_list)
# Specify dependencies for the library (if any)
target_link_libraries(plugin mints)

