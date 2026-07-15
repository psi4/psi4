include(TestingMacros)

if(ENABLE_Einsums OR ENABLE_Einsums2)
    add_regression_test(dfmp2f12-1 "psi;f12;df;quicktests")
endif()
