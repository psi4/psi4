include(TestingMacros)

add_regression_test(scf-cholesky-basis "psi;quicktests;scf;noc1")

# MKL_CBWR=AVX needed on macOS-Intel GHA runners (macos-14 upgrade changed FP results slightly)
if(APPLE)
    get_test_property(scf-cholesky-basis ENVIRONMENT _env)
    set_tests_properties(scf-cholesky-basis PROPERTIES ENVIRONMENT "${_env};MKL_CBWR=COMPATIBLE")
endif()
