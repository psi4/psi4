include(TestingMacros)

add_reference(nh3.npz ${CMAKE_CURRENT_BINARY_DIR})
add_regression_test(pcm_restart_pted "psi4;quicktests;pcm;ccsd")
