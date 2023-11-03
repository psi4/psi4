list(APPEND sources
  ein_test.cc
  )

if(TARGET Einsums::einsums)
    psi4_add_module(bin dummy_einsums sources)
    target_link_libraries(
      dummy_einsums
      PRIVATE
        Einsums::einsums
      )
endif()
