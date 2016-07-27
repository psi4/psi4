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
set(PLUGIN_LDFLAGS ${CMAKE_CXX_LINK_FLAGS} ${LIBC_INTERJECT})
add_definitions(-DPSI_TOP_SRCDIR=${PROJECT_ROOT}
                -Drestrict=${RESTRICT_KEYWORD}
)
set_source_files_properties(create_new_plugin.cc PROPERTIES
                -DPLUGIN_CXX=${CMAKE_CXX_COMPILER}
                -DPLUGIN_INCLUDES=${PSI4_ROOT}
                -DPLUGIN_FLAGS=${CMAKE_CXX_FLAGS}
                -DPLUGIN_DEFINES=""
                -DPLUGIN_OBJDIR=${CMAKE_BINARY_DIR}
                -DPLUGIN_LDFLAGS=${CMAKE_CXX_LINK_FLAGS}
)

#This library is special because we want this library to be a .so so that we
#can use the Python interface, even if the rest of Psi4 is static, but need to
#be able to build it statically for the binary
if(ENABLE_PYTHON_INTERFACE)
   add_library(psi4_main SHARED ${sources_list} $<TARGET_OBJECTS:versioned_code>)
 else()
   add_library(psi4_main STATIC ${sources_list} $<TARGET_OBJECTS:versioned_code>)
endif()

set(PSI4_EXTERNAL ${CMAKE_INSTALL_PREFIX}/external)
set(PSI4_RPATH ${CMAKE_INSTALL_PREFIX}/lib 
               ${PSI4_EXTERNAL}/boost/lib
               ${PSI4_EXTERNAL}/libderiv/lib
               ${PSI4_EXTERNAL}/libint/lib
               ${PSI4_EXTERNAL}/libefp/lib
)
set_target_properties(psi4_main PROPERTIES INSTALL_RPATH "${PSI4_RPATH}"
                                           PREFIX ""
                                           OUTPUT_NAME psi4
)                              

install(TARGETS psi4_main LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}
                          RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}
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
    ${LAPACK_LIBRARIES} ${PYTHON_LIBRARIES} ${Boost_LIBRARIES}
)
target_include_directories(psi4_main PUBLIC ${Boost_INCLUDE_DIRS} 
                                            ${LIBDERIV_INCLUDE_DIRS}
)
