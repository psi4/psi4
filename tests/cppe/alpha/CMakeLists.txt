include(TestingMacros)
file(COPY pna_6w.pot DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/../cppe-alpha)
add_regression_test(cppe-alpha "psi;cppe;addon;scf")

