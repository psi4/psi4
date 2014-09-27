set(SRC cleanup.cc f_act.cc hessian.cc params.cc ref_orbs.cc step.cc detcas.cc get_mo_info.cc indpairs.cc read_dens.cc rotate_orbs.cc test_bfgs.cc diis.cc gradient.cc ints.cc read_lag.cc setup_io.cc thetas.cc lagcalc.cc)
add_library(detcas ${SRC})
add_dependencies(detcas mints)
