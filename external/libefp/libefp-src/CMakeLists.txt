cmake_minimum_required(VERSION 3.0)
project(libefp LANGUAGES C)
set(sources_folder ${CMAKE_CURRENT_SOURCE_DIR}/src)
# List of sources
set(sources_list ${sources_folder}/ai_disp.c ${sources_folder}/balance.c 
        ${sources_folder}/bvec.c    ${sources_folder}/clapack.c 
        ${sources_folder}/disp.c    ${sources_folder}/efp.c 
        ${sources_folder}/elec.c    ${sources_folder}/elec_terms.c 
        ${sources_folder}/int.c     ${sources_folder}/log.c 
        ${sources_folder}/parse.c   ${sources_folder}/pol.c 
        ${sources_folder}/pol_direct.c ${sources_folder}/stream.c 
        ${sources_folder}/swf.c     ${sources_folder}/util.c 
        ${sources_folder}/xr.c
)

# Build static library
add_library(efp STATIC ${sources_list})
set_target_properties(efp PROPERTIES POSITION_INDEPENDENT_CODE ${BUILD_FPIC}
                                     COMPILE_FLAGS "-O1 -std=c99")

configure_file("libefpConfig.cmake.in" "libefpConfig.cmake" @ONLY)
install(FILES ${CMAKE_BINARY_DIR}/libefpConfig.cmake 
        DESTINATION ${CMAKE_INSTALL_PREFIX}
)
install(TARGETS efp DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/src/efp.h 
        DESTINATION ${CMAKE_INSTALL_PREFIX}/include/efp
)

