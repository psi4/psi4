include(TestingMacros)
file(COPY grid.dat DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
add_regression_test(scf-property "psi;quicktests;smoketests;scf;properties")
