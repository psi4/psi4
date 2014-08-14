set(SRC amp_write.cc c_clean.cc cache.cc cc2_hbar_extra.cc cc2_sigma.cc 
        cc3_HC1.cc cc3_HC1ET1.cc cceom.cc check_sum.cc dgeev_eom.cc 
        diag.cc diagSS.cc FDD.cc follow_root.cc form_diagonal.cc FSD.cc
        get_eom_params.cc get_moinfo.cc get_params.cc hbar_extra.cc 
        hbar_norms.cc local.cc local_guess.cc norm.cc norm_HC1.cc
        precondition.cc restart.cc restart_with_root.cc rzero.cc schmidt_add.cc
        sigma_full.cc sigmaCC3.cc sigmaCC3_RHF.cc sigmaDD.cc sigmaDS.cc
        sigmaSD.cc sigmaSS.cc sort_amps.cc sort_C.cc WabefDD.cc WabejDS.cc
        WamefSD.cc WbmfeDS.cc WmaijDS.cc WmbejDD.cc WmnefDD.cc 
        WmnieSD.cc WmnijDD.cc WnmjeDS.cc write_Rs.cc overlap.cc)
add_library(cceom ${SRC})
add_dependencies(cceom mints)
