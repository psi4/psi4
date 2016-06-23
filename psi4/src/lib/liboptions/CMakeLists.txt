set(sources_list print.cc 
                 python.cc 
                 liboptions.cc 
)
add_library(options ${sources_list})
target_link_libraries(options psi4util ${Boost_LIBRARIES})
target_include_directories(options PUBLIC ${Boost_INCLUDE_DIRS})
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} 
        DESTINATION psi4/src/lib/liboptions
        FILES_MATCHING PATTERN "*.h")
