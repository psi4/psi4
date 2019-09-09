include(TestingMacros)
file(COPY pna_6w.pot DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
add_regression_test(pe-scf-hf "psi;cppe;addon;scf")