set(sources_list #dealias.cc
                 dftensor.cc
                 #pstensor.cc
                 denominator.cc
                 pseudotrial.cc
                 fitter.cc
                 #schwarz.cc
                 fittingmetric.cc
                 #qr.cc
                 cholesky.cc
)
psi4_add_module(lib 3index sources_list mints fock)
