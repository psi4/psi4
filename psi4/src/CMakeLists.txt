### >> Go into psi4 subdirectory to compile libraries and modules <<
add_subdirectory(psi4)

####################################################################
### >> Merge libraries and modules into master psi4 library <<
### >> We also add in versioning, settings for plugins, and
### >> exports for Python here.

# Run the git version script, and update gitversion.h, if it's changed
add_custom_target(update_version
        COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/gitversion.py ${PROJECT_SOURCE_DIR}
        COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_CURRENT_BINARY_DIR}/gitversion.h.tmp ${CMAKE_CURRENT_BINARY_DIR}/gitversion.h
        COMMAND ${CMAKE_COMMAND} -E remove ${CMAKE_CURRENT_BINARY_DIR}/gitversion.h.tmp
        COMMENT "Generating Git info")

add_library(versioned_code OBJECT version.cc python.cc psi_start.cc)
add_dependencies(versioned_code update_version)
set_target_properties(versioned_code PROPERTIES POSITION_INDEPENDENT_CODE 1)
target_include_directories(versioned_code PUBLIC ${PYBIND11_INCLUDE_DIRS} ${CMAKE_CURRENT_BINARY_DIR})

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

# TODO: Having PSI4_ROOT for plugins isn't right
set(PLUGIN_INCLUDES "-I${PSI4_ROOT} -I${PYTHON_INCLUDE_DIR} ")

add_definitions(-DPSI_TOP_SRCDIR=${PROJECT_ROOT}
        -Drestrict=${RESTRICT_KEYWORD}
        -DPLUGIN_CXX=${CMAKE_CXX_COMPILER}
        -DPLUGIN_INCLUDES=${PLUGIN_INCLUDES}
        -DPLUGIN_FLAGS="${CMAKE_CXX_FLAGS}"
        -DPLUGIN_DEFINES="-fPIC"
        -DPLUGIN_OBJDIR=${CMAKE_BINARY_DIR}
        -DPLUGIN_LDFLAGS="${CMAKE_CXX_LINK_FLAGS}"
        )

set(psi4_all_sources ${sources_list} $<TARGET_OBJECTS:versioned_code>)
add_library(psi4 STATIC ${psi4_all_sources})
if(TARGET gdma::gdma)
    target_compile_definitions(versioned_code PRIVATE $<TARGET_PROPERTY:gdma::gdma,INTERFACE_COMPILE_DEFINITIONS>)
endif()
if(TARGET CheMPS2::chemps2)
    target_compile_definitions(versioned_code PRIVATE $<TARGET_PROPERTY:CheMPS2::chemps2,INTERFACE_COMPILE_DEFINITIONS>)
endif()

install(TARGETS psi4 ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR})

#set_property(GLOBAL APPEND PROPERTY LIBLIST psi4)
#target_link_libraries(psi4 PRIVATE mints cubeprop plugin parallel detci
#        ccresponse ccdensity cchbar fnocc cctriples ccenergy cclambda cceom
#        cctransort ccsort fisapt dfmp2 dcft findif psimrcc dfocc efp_interface adc
#        mcscf occ optking functional scf_solver scfgrad transqt2 moinfo 3index
#        thermo mrcc scf sapt gdma_interface dmrg
#        )
get_property(psi4_libraries GLOBAL PROPERTY LIBLIST)

target_link_libraries(psi4 PRIVATE ${psi4_libraries})
target_link_libraries(psi4 PRIVATE ${PYTHON_LIBRARIES})
target_link_libraries(psi4 PRIVATE ${CMAKE_THREAD_LIBS_INIT})
target_link_libraries(psi4 PRIVATE ${LAPACK_LIBRARIES})
target_link_libraries(psi4 PRIVATE ${BLAS_LIBRARIES})

target_include_directories(psi4 PUBLIC ${PYBIND11_INCLUDE_DIRS})
target_include_directories(psi4 PRIVATE $<TARGET_PROPERTY:libint::libint,INTERFACE_INCLUDE_DIRECTORIES>)  # comp def instead? only b/c boost export here not w/i mints


### >>  Stuff below this is for RPaths and doesn't work quite right

# TODO: Investiage relative RPaths.
# Add both the build and install RPaths so we can test from build dir
set(PSI4_EXTERNAL ${CMAKE_INSTALL_PREFIX}/external)
set(PSI4_BIN_EXTERNAL ${CMAKE_BINARY_DIR}/stage/${PSI4_EXTERNAL})
foreach(ex_proj libderiv libint libefp)
    set(PSI4_RPATH
            "${PSI4_RPATH};${PSI4_EXTERNAL}/${ex_proj}/lib;${PSI4_BIN_EXTERNAL}/${ex_proj}/lib")
endforeach()

# The other big source of RPaths are the math libraries
foreach(math_comp ${LAPACK_LIBRARIES};${BLAS_LIBRARIES})
    get_filename_component(math_lib ${math_comp} DIRECTORY)
    set(PSI4_RPATH "${PSI4_RPATH};${math_lib}")
endforeach()

add_executable(psi4-bin psi4.cc)
install(TARGETS psi4-bin RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_BINDIR})

get_filename_component(PYTHON_LIB_DIR ${PYTHON_LIBRARIES} DIRECTORY)

set(psi4_install_lib ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR})
set(PSI4_RPATH ${psi4_install_lib}
        ${DESTDIR}/${psi4_install_lib}
        ${PYTHON_LIB_DIR}
        )
set_target_properties(psi4-bin PROPERTIES OUTPUT_NAME psi4
        INSTALL_RPATH "${PSI4_RPATH}"
        INSTALL_RPATH_USE_LINK_PATH TRUE
        BUILD_WITH_INSTALL_RPATH TRUE
        )

message("LIBS: PYTHON   ${PYTHON_LIBRARIES}")
message("LIBS: PSI4     ${psi4_libraries}")
message("LIBS: CHEMPS2  ${CHEMPS2_LIBRARIES}")
message("LIBS: GDMA")
message("LIBS: DERIV    ${LIBDERIV_LIBRARIES}")
message("LIBS: THREAD   ${CMAKE_THREAD_LIBS_INIT}")
message("LIBS: LAPACK   ${LAPACK_LIBRARIES}")
message("LIBS: BLAS     ${BLAS_LIBRARIES}")
message("LIBS: Fortran  ${CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES}")

target_link_libraries(psi4-bin PRIVATE psi4)

