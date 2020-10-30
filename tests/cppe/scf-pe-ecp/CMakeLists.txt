include(TestingMacros)

file(COPY potential.pot DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/../cppe-scf-pe-ecp)
file(COPY potential_sulfur.pot DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/../cppe-scf-pe-ecp)
add_regression_test(cppe-scf-pe-ecp "psi;cppe;addon;scf")
