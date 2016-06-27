set(sources_list print.cc 
                 python.cc 
                 liboptions.cc 
)
psi4_add_library(options sources_list)
target_link_libraries(options psi4util ${Boost_LIBRARIES})
target_include_directories(options PUBLIC ${Boost_INCLUDE_DIRS})
