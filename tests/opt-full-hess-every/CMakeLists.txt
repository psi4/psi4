include(TestingMacros)

if(MAX_AM_ERI GREATER_EQUAL 6)
    add_regression_test(opt-full-hess-every "psi;quicktests;opt")
endif()
