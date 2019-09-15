include(TestingMacros)

file(COPY pna_6w.pot DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/../cppe-scf-hf)
add_regression_test(cppe-scf-hf "psi;cppe;addon;scf")
