set(SRC amps.cc build_A.cc build_X.cc cache.cc check_energy.cc deanti.cc dipole.cc energy.cc fold.cc get_moinfo.cc get_params.cc Iab.cc Iai.cc Iia.cc Iij.cc lag.cc mp2.cc opdm.cc priority.cc relax_I.cc relax_opdm.cc sort_amps.cc sort_I.cc sort_opdm.cc sort_twopdm.cc twopdm.cc write_data.cc Zvector.cc)
add_library(mp2 ${SRC})
add_dependencies(mp2 mints)
