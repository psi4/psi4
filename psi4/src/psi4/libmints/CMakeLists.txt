set(sources_list local.cc 
                 onebody.cc 
                 x2cint.cc 
                 orbitalspace.cc 
                 osrecur.cc 
                 maketab.cc 
                 efpmultipolepotential.cc 
                 rel_potential.cc 
                 oeprop.cc 
                 writer.cc 
                 transform.cc 
                 sieve.cc 
                 multipolesymmetry.cc 
                 shellrotation.cc 
                 deriv.cc 
                 overlap.cc 
                 integralparameters.cc 
                 twobody.cc 
                 vector.cc 
                 sobasis.cc 
                 view.cc 
                 cartesianiter.cc 
                 basisset.cc 
                 electrostatic.cc 
                 wavefunction.cc 
                 basisset_parser.cc 
                 irrep.cc 
                 eribase.cc 
                 fjt.cc 
                 potentialint.cc 
                 chartab.cc 
                 corrtab.cc 
                 quadrupole.cc 
                 eri.cc 
                 symop.cc 
                 benchmark.cc 
                 get_writer_file_prefix.cc 
                 3coverlap.cc 
                 petitelist.cc 
                 solidharmonics.cc 
                 orthog.cc 
                 electricfield.cc 
                 multipoles.cc 
                 dipole.cc 
                 sointegral.cc 
                 extern.cc 
                 nabla.cc 
                 factory.cc 
                 psimath.cc 
                 dimension.cc 
                 molecule.cc 
                 intvector.cc 
                 potential.cc 
                 mintshelper.cc 
                 coordentry.cc 
                 kinetic.cc 
                 tracelessquadrupole.cc 
                 pseudospectral.cc 
                 integral.cc 
                 matrix.cc 
                 svd.cc 
                 gshell.cc 
                 integraliter.cc 
                 pointgrp.cc 
                 rep.cc 
                 cdsalclist.cc 
                 erd_eri.cc 
                 angularmomentum.cc
)
add_definitions("-Drestrict=${RESTRICT_KEYWORD}")
psi4_add_library(mints sources_list iwl options psi4util trans)
target_link_libraries(mints PRIVATE libint::libint)
if(TARGET dkh::dkh)
    target_link_libraries(mints PRIVATE dkh::dkh)
endif()
if(TARGET liberd::erd)
    target_link_libraries(mints PRIVATE liberd::erd)
endif()

