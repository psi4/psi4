cmake_minimum_required(VERSION 3.0)
project(libefp LANGUAGES C)
set(sources_list "")
# List of sources
set(tmp ai_disp.c 
         balance.c 
         bvec.c 
         clapack.c 
         disp.c 
         efp.c 
         elec.c 
         elec_terms.c 
         int.c 
         log.c 
         parse.c 
         pol.c 
         pol_direct.c 
         stream.c 
         swf.c 
         util.c 
         xr.c
)

foreach(fl ${tmp})
    list(APPEND sources_list ${CMAKE_CURRENT_SOURCE_DIR}/src/${fl})
endforeach(fl)

# Build static library
add_library(efp STATIC ${sources_list})
set_property(TARGET efp PROPERTY COMPILE_FLAGS "-O1 -std=c99")

configure_file("libefpConfig.cmake.in" "libefpConfig.cmake" @ONLY)
install(FILES ${CMAKE_BINARY_DIR}/libefpConfig.cmake 
        DESTINATION ${CMAKE_INSTALL_PREFIX}
)
install(TARGETS efp DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/src/efp.h 
        DESTINATION ${CMAKE_INSTALL_PREFIX}/include/efp
)

