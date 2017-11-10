include(TestingMacros)

add_reference(nh3.npz ${CMAKE_CURRENT_BINARY_DIR})
add_regression_test(pcmsolver-ccsd-pted-restart "psi;quicktests;pcmsolver;addon;ccsd;pted")
