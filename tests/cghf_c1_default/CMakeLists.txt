include(TestingMacros)

if(ENABLE_Einsums)
    add_regression_test(cghf_c1_default "psi;quicktests;smoketests;cghf")
endif()
