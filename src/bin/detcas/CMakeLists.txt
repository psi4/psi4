set(SRC cleanup.cc detcas.cc diis.cc f_act.cc get_mo_info.cc gradient.cc
hessian.cc indpairs.cc ints.cc lagcalc.cc params.cc read_dens.cc read_lag.cc
ref_orbs.cc rotate_orbs.cc setup_io.cc step.cc test_bfgs.cc thetas.cc)
add_library(detcas ${SRC})
add_dependencies(detcas mints)
