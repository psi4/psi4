include(TestingMacros)

file(COPY ../scf-hf/pna_6w.pot DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/../cppe-tdscf)
add_regression_test(cppe-tdscf "psi;cppe;addon;tdscf")
