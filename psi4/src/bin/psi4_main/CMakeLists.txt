# Run the git version script, and update gitversion.h, if it's changed
add_custom_target(update_version
    COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/gitversion.py ${PROJECT_SOURCE_DIR}
    COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_BINARY_DIR}/gitversion.h.tmp ${CMAKE_CURRENT_BINARY_DIR}/gitversion.h
    COMMAND ${CMAKE_COMMAND} -E remove ${CMAKE_CURRENT_BINARY_DIR}/gitversion.h.tmp
    COMMENT "Generating Git info")

add_library(versioned_code OBJECT version.cc python.cc psi_start.cc)
add_dependencies(versioned_code update_version)
set_target_properties(versioned_code PROPERTIES POSITION_INDEPENDENT_CODE 1)
target_include_directories(versioned_code PUBLIC ${Boost_INCLUDE_DIRS} ${CMAKE_CURRENT_BINARY_DIR})
 
set(sources_list export_psio.cc 
                 export_mints.cc 
                 psi_stop.cc 
                 export_functional.cc 
                 export_oeprop.cc 
                 export_plugins.cc 
                 export_blas_lapack.cc 
                 export_benchmarks.cc 
                 export_efp.cc 
                 export_cubeprop.cc 
                 clean.cc 
                 create_new_plugin.cc 
                 script.cc 
                 set_memory.cc 
                 read_options.cc 
)

set(PLUGIN_INCLUDES "-I${PSI4_ROOT} -I${Boost_INCLUDE_DIR} -I${PYTHON_INCLUDE_DIR} ")

add_definitions(-DPSI_TOP_SRCDIR=${PROJECT_ROOT}
    -Drestrict=${RESTRICT_KEYWORD}
    -DPLUGIN_CXX=${CMAKE_CXX_COMPILER}
    -DPLUGIN_INCLUDES=${PLUGIN_INCLUDES}
    -DPLUGIN_FLAGS="${CMAKE_CXX_FLAGS}"
    -DPLUGIN_DEFINES="-fPIC"
    -DPLUGIN_OBJDIR=${CMAKE_BINARY_DIR}
    -DPLUGIN_LDFLAGS="${CMAKE_CXX_LINK_FLAGS}"
)

#This library is special because we want this library to be a .so so that we
#can use the Python interface, even if the rest of Psi4 is static, but need to
#be able to build it statically for the binary
set(psi4_all_sources ${sources_list} $<TARGET_OBJECTS:versioned_code>)
if(ENABLE_PYTHON_INTERFACE)
   add_library(psi4_main SHARED ${psi4_all_sources})
 else()
   add_library(psi4_main STATIC ${psi4_all_sources})
endif()

set(psi4_destination ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR})

install(TARGETS psi4_main LIBRARY DESTINATION ${psi4_destination}
                          RUNTIME DESTINATION ${psi4_destination}
)

install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} 
        DESTINATION psi4/src/bin/psi4_main
        FILES_MATCHING PATTERN "*.h" PATTERN "*.hpp"
)

set_property(GLOBAL APPEND PROPERTY LIBLIST psi4_main)
target_link_libraries(psi4_main PRIVATE  mints cubeprop plugin parallel detci
    ccresponse ccdensity cchbar fnocc cctriples ccenergy cclambda cceom
    cctransort ccsort fisapt dfmp2 dcft findif psimrcc dfocc efp_interface adc
    mcscf occ optking functional scf_solver scfgrad transqt2 moinfo 3index deriv 
    int thermo mrcc scf sapt gdma_interface dmrg
    ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} ${PYTHON_LIBRARIES} ${Boost_LIBRARIES}
)
target_include_directories(psi4_main PUBLIC ${Boost_INCLUDE_DIRS} 
                                            ${LIBDERIV_INCLUDE_DIRS}
)


################################################################################
#  Stuff below this is for RPaths and doesn't work quite right
################################################################################


if(BUILD_SHARED_LIBS)
   #Add both the build and install RPaths so we can test from build dir
   set(PSI4_EXTERNAL ${CMAKE_INSTALL_PREFIX}/external)
   set(PSI4_BIN_EXTERNAL ${CMAKE_BINARY_DIR}/stage/${PSI4_EXTERNAL})
   foreach(ex_proj boost libderiv libint libefp)
      set(PSI4_RPATH 
         "${PSI4_RPATH};${PSI4_EXTERNAL}/${ex_proj}/lib;${PSI4_BIN_EXTERNAL}/${ex_proj}/lib")
   endforeach()
endif()

#The other big source of RPaths are the math libraries
foreach(math_comp "${LAPACK_LIBRARIES} ${BLAS_LIBRARIES}")
    get_filename_component(math_lib ${math_comp} DIRECTORY)
    set(PSI4_RPATH "${PSI4_RPATH};${math_lib}")
endforeach()


set_target_properties(psi4_main PROPERTIES INSTALL_RPATH "${PSI4_RPATH}"
                                           INSTALL_RPATH_USE_LINK_PATH TRUE
                                           BUILD_WITH_INSTALL_RPATH TRUE
                                           PREFIX ""
                                           OUTPUT_NAME psi4
)                              


