set(SRC amp_write.cc analyze.cc cache.cc cc2_hbar_extra.cc cc2_LHX1Y1.cc cc2_LHX1Y2.cc cc2_sort_X.cc cc2_X1.cc cc2_X2.cc ccresponse.cc compute_X.cc converged.cc denom.cc diis.cc get_moinfo.cc get_params.cc hbar_extra.cc HXY.cc init_X.cc lambda_residuals.cc LCX.cc LHX1Y1.cc LHX1Y2.cc LHX2Y2.cc linresp.cc local.cc optrot.cc pertbar.cc polar.cc preppert.cc print_X.cc pseudopolar.cc roa.cc save_X.cc sort_lamps.cc sort_pert.cc sort_X.cc update_X.cc X1.cc X2.cc)
add_library(ccresponse ${SRC})
add_dependencies(ccresponse mints)
