set(sources_list quad.cc 
                 stl_string.cc 
                 memory_manager.cc 
                 exception.cc 
                 combinations.cc 
)

# Build static library
add_library(psi4util ${sources_list})
target_link_libraries(psi4util options)
include_directories(${Boost_INCLUDE_DIRS})
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} 
        DESTINATION psi4/src/lib/libpsi4util
        FILES_MATCHING PATTERN "*.h")
