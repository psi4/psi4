find_package(libderiv REQUIRED COMPONENTS ${LIBINT_OPT_AM})

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

if(ENABLE_DKH)
   list(APPEND sources_list dkh2-dkh4_main.F90)
endif()
add_library(mints ${sources_list})
target_include_directories(mints PUBLIC deriv ${Boost_INCLUDE_DIRS})
target_link_libraries(mints qt iwl deriv ${Boost_LIBRARIES})

install(TARGETS mints DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} DESTINATION psi4/src/lib/libmints
        FILES_MATCHING PATTERN "*.h")
