include(TestingMacros)

add_regression_test(dft-water-dimer "psi;dft;scf")
set_tests_properties(dft-water-dimer PROPERTIES TIMEOUT 3000) 
