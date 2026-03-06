list(APPEND sources
  fisapt.cc
  local2.cc
  )
psi4_add_module(bin fisapt sources)

if(TARGET gauxc::gauxc)
  target_include_directories(fisapt
    PUBLIC
      $<TARGET_PROPERTY:gauxc::gauxc,INTERFACE_INCLUDE_DIRECTORIES>
    )
endif()


# NOTE: logic taken from f12 CMakeLists.txt
if(TARGET Einsums::Einsums)
    if(DEFINED ENV{CONDA_TOOLCHAIN_BUILD})
        if(APPLE)
            target_compile_definitions(
              fisapt
              PUBLIC
                _LIBCPP_DISABLE_AVAILABILITY
                H5CPP_USE_OMP_ALIGNED_ALLOC
              )
        endif()
    endif()
    target_link_libraries(fisapt PUBLIC diis Einsums::Einsums)
else()
    target_link_libraries(fisapt PUBLIC diis)
endif()

# target_link_libraries(fisapt PUBLIC diis)
