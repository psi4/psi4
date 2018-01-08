include(TestingMacros)

file(COPY Ne.6311G.INTDUMP.ref DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
add_regression_test(python-fcidump "psi;quicktests;smoketests;fcidump;python")
