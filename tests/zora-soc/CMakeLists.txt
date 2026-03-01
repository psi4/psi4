include(TestingMacros)

if(ENABLE_Einsums)
    add_regression_test(zora-soc "psi;quicktests;zora")
endif()
