### >> Go into psi4 subdirectory to compile libraries and modules <<
add_subdirectory(psi4)

####################################################################
### >> Merge libraries and modules into master psi4 library <<
### >> We also add in versioning, settings for plugins, and
### >> exports for Python here.

set(sources_list export_psio.cc
        export_mints.cc
        export_fock.cc
        export_functional.cc
        export_oeprop.cc
        export_plugins.cc
        export_blas_lapack.cc
        export_benchmarks.cc
        export_efp.cc
        export_cubeprop.cc
        export_misc.cc
        create_new_plugin.cc
        read_options.cc
        core.cc
        )

add_library(core SHARED ${sources_list})

get_property(psi4_binmodules GLOBAL PROPERTY BINLIST)
target_link_libraries(core PRIVATE ${PRE_LIBRARY_OPTION} dpd qt ${POST_LIBRARY_OPTION})

if(TARGET dkh::dkh)
    target_compile_definitions(core PRIVATE $<TARGET_PROPERTY:dkh::dkh,INTERFACE_COMPILE_DEFINITIONS>)
endif()
target_link_libraries(core PRIVATE ${psi4_binmodules})
target_link_libraries(core PRIVATE ${LIBC_INTERJECT})
target_link_libraries(core PUBLIC pybind11::module)
target_link_libraries(core PRIVATE Threads::Threads)
# LAPACK & BLAS linking attached to modules in BIN/LIBLIST to maximally defer

target_include_directories(core PRIVATE $<TARGET_PROPERTY:libint::libint,INTERFACE_INCLUDE_DIRECTORIES>)  # comp def instead? only b/c boost export here not w/i mints
target_include_directories(core INTERFACE $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

install(TARGETS core
        EXPORT "psi4Targets"
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}${PYMOD_INSTALL_LIBDIR}/psi4
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
install(EXPORT "psi4Targets"
        NAMESPACE "psi4::"
        DESTINATION ${CMAKECONFIG_INSTALL_DIR})

### >>  Stuff below this is for RPaths and doesn't work quite right

# TODO: Investiage relative RPaths.
# Add both the build and install RPaths so we can test from build dir
#set(PSI4_EXTERNAL ${CMAKE_INSTALL_PREFIX}/external)
#set(PSI4_BIN_EXTERNAL ${CMAKE_BINARY_DIR}/stage${PSI4_EXTERNAL})
#foreach(ex_proj libderiv libint libefp)
#    set(PSI4_RPATH
#            "${PSI4_RPATH};${PSI4_EXTERNAL}/${ex_proj}/lib;${PSI4_BIN_EXTERNAL}/${ex_proj}/lib")
#endforeach()
set(PSI4_RPATH "${CMAKE_BINARY_DIR}/stage${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}")

# The other big source of RPaths are the math libraries
#foreach(math_comp tgt::lapack)
#    message("${math_comp}")
#    get_filename_component(math_lib ${math_comp} DIRECTORY)
#    set(PSI4_RPATH "${PSI4_RPATH};${math_lib}")
#endforeach()

get_filename_component(PYTHON_LIB_DIR ${PYTHON_LIBRARIES} DIRECTORY)

set(psi4_install_lib ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR})
set(PSI4_RPATH ${psi4_install_lib}
        ${DESTDIR}/${psi4_install_lib}
        ${PYTHON_LIB_DIR}
        )

set_target_properties(core PROPERTIES PREFIX "${PYTHON_MODULE_PREFIX}"  # for python module
                           OUTPUT_NAME core
                           EXPORT_NAME core
                           SUFFIX ".so"  # for python module
                           #SUFFIX "${PYTHON_MODULE_EXTENSION}"  # for python module
                           INSTALL_RPATH "${PSI4_RPATH}"
                           INSTALL_RPATH_USE_LINK_PATH TRUE
                           BUILD_WITH_INSTALL_RPATH TRUE)

