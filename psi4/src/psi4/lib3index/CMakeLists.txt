set(sources_list dftensor.cc
                 dfhelper.cc
                 denominator.cc
                 fittingmetric.cc
                 cholesky.cc
)
psi4_add_module(lib 3index sources_list mints fock)
