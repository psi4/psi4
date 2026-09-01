include(TestingMacros)

if(ENABLE_Einsums)
    add_regression_test(cghf_without_diis "psi;quicktests;cghf")
endif()
