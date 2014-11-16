set(headers_list "")
# List of headers
list(APPEND headers_list electrostatic.h writer_file_prefix.h gridblock.h factory.h wavefunction.h oeprop.h cubefile.h potentialint.h benchmark.h overlap.h cdsalclist.h angularmomentum.h vector3.h potential.h local.h integral.h sointegral_onebody.h dimension.h 3coverlap.h corrtab.h vector.h extern.h pointgrp.h molecule.h sieve.h sointegral.h quadrupole.h typedefs.h electricfield.h fjt.h psimath.h osrecur.h sobasis.h integralparameters.h sointegral_twobody.h mints.h onebody.h petitelist.h efpmultipolepotential.h coordentry.h orbitalspace.h nabla.h eri.h orthog.h serializers.h erd_eri.h tracelessquadrupole.h basisset_parser.h dipole.h rel_potential.h kinetic.h dcd.h twobody.h basisset.h shellrotation.h gshell.h view.h multipoles.h mintshelper.h writer.h pybuffer.h cartesianiter.h multipolesymmetry.h deriv.h pseudospectral.h matrix.h )

# If you want to remove some headers specify them explictly here
if(DEVELOPMENT_CODE)
   list(REMOVE_ITEM headers_list "")
else()
   list(REMOVE_ITEM headers_list "")
endif()
# Sort alphabetically
list(SORT headers_list)

set(sources_list "")
# List of sources
list(APPEND sources_list local.cc onebody.cc orbitalspace.cc osrecur.cc maketab.cc efpmultipolepotential.cc rel_potential.cc oeprop.cc writer.cc transform.cc cubefile.cc sieve.cc multipolesymmetry.cc shellrotation.cc deriv.cc overlap.cc integralparameters.cc twobody.cc vector.cc sobasis.cc view.cc cartesianiter.cc basisset.cc electrostatic.cc wavefunction.cc basisset_parser.cc irrep.cc eribase.cc fjt.cc potentialint.cc chartab.cc corrtab.cc quadrupole.cc gridprop.cc eri.cc symop.cc benchmark.cc get_writer_file_prefix.cc 3coverlap.cc petitelist.cc solidharmonics.cc orthog.cc electricfield.cc multipoles.cc dipole.cc sointegral.cc extern.cc nabla.cc factory.cc psimath.cc dimension.cc molecule.cc intvector.cc potential.cc mintshelper.cc coordentry.cc kinetic.cc tracelessquadrupole.cc pseudospectral.cc integral.cc matrix.cc svd.cc gshell.cc integraliter.cc pointgrp.cc rep.cc cdsalclist.cc erd_eri.cc angularmomentum.cc )

# If you want to remove some sources specify them explictly here
if(DEVELOPMENT_CODE)
   list(REMOVE_ITEM sources_list "gridprop.cc")
else()
   list(REMOVE_ITEM sources_list "gridprop.cc")
endif()

# Write list of files to be passed to cloc for counting lines of code.
# Only files that are actually compiled are counted.
set(to_count "${sources_list}" "${headers_list}")
write_to_cloc_list("${to_count}")

# Build static library
add_library(mints STATIC ${sources_list})
# Specify dependencies for the library (if any)
add_dependencies(mints int deriv)
set(libs_to_merge mints ${libs_to_merge} PARENT_SCOPE)
if(BUILD_CUSTOM_BOOST)
   add_dependencies(mints custom_boost)
endif()
install(TARGETS mints ARCHIVE DESTINATION lib)

# Sets install directory for all the headers in the list
install_list_FILES("${headers_list}" include/libmints)
