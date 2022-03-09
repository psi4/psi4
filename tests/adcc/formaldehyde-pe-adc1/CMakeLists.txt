include(TestingMacros)

file(COPY fa_6w.pot DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/../adcc-formaldehyde-pe-adc1)
add_regression_test(adcc-formaldehyde-pe-adc1 "psi;adc;addon;cppe;quicktests")
