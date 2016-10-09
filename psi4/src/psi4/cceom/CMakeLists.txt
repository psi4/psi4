set(sources_list local.cc form_diagonal.cc WmnieSD.cc WabejDS.cc diagSS.cc check_sum.cc sigma_full.cc hbar_extra.cc sigmaSS.cc restart_with_root.cc cache.cc sigmaCC3.cc rzero.cc overlap.cc sigmaCC3_RHF.cc local_guess.cc get_params.cc 
#read_guess.cc 
dgeev_eom.cc schmidt_add.cc norm_HC1.cc restart.cc sort_amps.cc cc3_HC1ET1.cc FSD.cc write_Rs.cc follow_root.cc WmnefDD.cc WamefSD.cc get_eom_params.cc sort_C.cc FDD.cc WmaijDS.cc precondition.cc WabefDD.cc get_moinfo.cc sigmaDD.cc WmbejDD.cc cc2_hbar_extra.cc cc3_HC1.cc c_clean.cc cceom.cc amp_write.cc cc2_sigma.cc WmnijDD.cc WbmfeDS.cc sigmaDS.cc sigmaSD.cc hbar_norms.cc norm.cc WnmjeDS.cc diag.cc )
psi4_add_module(bin cceom sources_list mints)
