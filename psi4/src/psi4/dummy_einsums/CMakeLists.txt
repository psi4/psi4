list(APPEND sources
  ein_test.cc
  )

if(TARGET Einsums::Einsums)
    psi4_add_module(bin dummy_einsums sources)
    target_link_libraries(
      dummy_einsums
      PRIVATE
        Einsums::Einsums
      )

    if(DEFINED ENV{CONDA_TOOLCHAIN_BUILD})
        if(APPLE)
            # * possibly these should be exported from Einsums somehow
            # * when using conda clang, compiling against an old SDK &
            #   an openmp symbols is absent
            # * afaik, envvar above is defined when `conda activate envvar_w_compilers`.
            #   There's a CONDA_BACKUP_CONDA_TOOLCHAIN_BUILD envvar when not activated.
            #   So hopefully, this doesn't interfere with system compilers.
            target_compile_definitions(
              dummy_einsums
              PUBLIC
                _LIBCPP_DISABLE_AVAILABILITY
                H5CPP_USE_OMP_ALIGNED_ALLOC
              )
        endif()
    endif()
endif()
