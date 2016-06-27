set(sources_list lreri.cc 
                 thcew.cc 
                 laplace.cc 
                 thce.cc 
)
psi4_add_library(thce sources_list)
target_link_libraries(thce mints)
