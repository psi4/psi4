include(TestingMacros)

if(ENABLE_Einsums)
    add_regression_test(cghf_spin_square "psi;quicktests;cghf")
endif()
