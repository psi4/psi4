set(SRC cache.cc count_ijk.cc ET_AAA.cc ET_AAB.cc ET_ABB.cc ET_BBB.cc ET_RHF.cc ET_UHF_AAA.cc ET_UHF_AAB.cc ET_UHF_ABB.cc ET_UHF_BBB.cc get_moinfo.cc T3_grad_RHF.cc T3_grad_UHF_AAA.cc T3_grad_UHF_AAB.cc T3_grad_UHF_BBA.cc T3_grad_UHF_BBB.cc T3_UHF_AAA.cc T3_UHF_AAB.cc T3_UHF_ABC.cc test_abc_loops.cc transpose_integrals.cc triples.cc)
add_library(cctriples ${SRC})
