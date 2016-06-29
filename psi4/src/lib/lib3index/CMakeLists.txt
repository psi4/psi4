set(sources_list dealias.cc dftensor.cc pstensor.cc denominator.cc pseudotrial.cc fitter.cc schwarz.cc fittingmetric.cc qr.cc cholesky.cc )
psi4_add_library(3index sources_list)
target_link_libraries(3index mints)
