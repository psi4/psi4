include(TestingMacros)
file(COPY grid.dat DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
add_regression_test(scf_property "psi;quicktests;scf;properties")
