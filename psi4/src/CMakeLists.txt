### >> Create core target, so we may link modules to it directly <<
### >> This allow to propagate dependencies upwards sanely
add_library(core SHARED "")

### >> Go into psi4 subdirectory to compile libraries and modules <<
add_subdirectory(psi4)

####################################################################
### >> Merge libraries and modules into master psi4 library <<
### >> We also add in versioning, settings for plugins, and
### >> exports for Python here.

target_sources(core
  PRIVATE
    export_benchmarks.cc
    export_blas_lapack.cc
    export_cubeprop.cc
    export_diis.cc
    export_fock.cc
    export_functional.cc
    export_mints.cc
    export_misc.cc
    export_oeprop.cc
    export_pcm.cc
    export_plugins.cc
    export_psio.cc
    export_trans.cc
    export_wavefunction.cc
    export_options.cc
    create_new_plugin.cc
    read_options.cc
    python_data_type.cc
    core.cc
  )

set_target_properties(core
  PROPERTIES
    CXX_VISIBILITY_PRESET hidden
    POSITION_INDEPENDENT_CODE 1
    VISIBILITY_INLINES_HIDDEN 1
  )

get_property(psi4_binmodules GLOBAL PROPERTY LINKBINLIST)
get_property(psi4_libmodules GLOBAL PROPERTY LINKLIBLIST)
target_link_libraries(core
  PRIVATE
    ${psi4_binmodules}
    ${PRE_LIBRARY_OPTION}
    ${psi4_libmodules}
    ${POST_LIBRARY_OPTION}
    ${LIBC_INTERJECT}
    tgt::lapack # OpenMP flags to core.cc
  PUBLIC
    pybind11::module
  )

if(Fortran_ENABLED AND CMAKE_Fortran_COMPILER_ID MATCHES Intel)
  # Enable call to for_rtl_init_() which is required if using the
  # Intel fortran compiler
  target_compile_definitions(core
    PRIVATE
      INTEL_Fortran_ENABLED
    )
endif()

if(MSVC)
  target_link_libraries(core PRIVATE Ws2_32)
  # Increase the number of addressable sections
  target_compile_options(core PRIVATE "/bigobj")
endif()

target_include_directories(core
  INTERFACE
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
  )

install(TARGETS core
        EXPORT "psi4Targets"
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}${PYMOD_INSTALL_LIBDIR}/psi4
        RUNTIME DESTINATION ${CMAKE_INSTALL_LIBDIR}${PYMOD_INSTALL_LIBDIR}/psi4)
install(EXPORT "psi4Targets"
        NAMESPACE "psi4::"
        DESTINATION ${CMAKECONFIG_INSTALL_DIR})

# <<<  RPATH  >>>

set(_full_path_staged_libdir "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}")
file(RELATIVE_PATH _rel ${_full_path_staged_libdir}${PYMOD_INSTALL_LIBDIR}/psi4
                        ${_full_path_staged_libdir})
if(APPLE)
    set(psi4_RPATH "@loader_path/${_rel}")
else()
    set(psi4_RPATH "$ORIGIN/${_rel}")
endif()

foreach(_al ${_addons})
    get_filename_component(_ald ${_al} DIRECTORY)

    if(${_al} MATCHES "-Wl,")
        list(REMOVE_ITEM _addons ${_al})
        #message("rpath: ${_al} dropped because linker directive")
        continue()
    endif()

    if(${_al} MATCHES "\.${CMAKE_STATIC_LIBRARY_SUFFIX}$")
        list(REMOVE_ITEM _addons ${_al})
        #message("rpath: ${_al} dropped because static")
        continue()
    endif()

    if(NOT ${_ald} STREQUAL "")
        list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES ${_ald} _in_cpild)
        if(_in_cpild GREATER -1)
            list(REMOVE_ITEM _addons ${_al})
            #message("rpath: ${_al} dropped because system lib")
            continue()
        endif()

        if(${_ald} STREQUAL ${_full_path_staged_libdir})
            list(REMOVE_ITEM _addons ${_al})
            #message("rpath: ${_al} dropped because internal relative to core.so")
            continue()
        endif()

        list(APPEND psi4_RPATH ${_ald})
        #message("rpath: ${_al} included")
    endif()
endforeach()

list(REMOVE_DUPLICATES psi4_RPATH)
message(STATUS "Psi4 rpath: ${psi4_RPATH}")

set_target_properties(core PROPERTIES PREFIX "${PYTHON_MODULE_PREFIX}"  # for python module
                           OUTPUT_NAME core
                           EXPORT_NAME core
                           SUFFIX "${PYTHON_MODULE_EXTENSION}" # for python module
                           INSTALL_RPATH "${psi4_RPATH}"
                           BUILD_WITH_INSTALL_RPATH ON)
