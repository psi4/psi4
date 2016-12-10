set(sources_list #zmat_point.cc
                 lapack_intfc.cc
                 #slaterdset.cc
                 #strncpy.cc
                 #sort.cc
                 #mat_in.cc
                 dx_write.cc
                 #orient_fragment.cc
                 #rotate_vecs.cc
                 dirprd_block.cc
                 cc_wfn.cc pople.cc
                 schmidt_add.cc
                 mat_print.cc
                 filter.cc
                 ras_set.cc
                 fill_sym_matrix.cc
                 #rootfind.cc
                 cc_excited.cc
                 probabil.cc
                 dot_block.cc
                 timer.cc
                 dx_read.cc
                 blas_intfc.cc
                 normalize.cc
                 newmm_rking.cc
                 3d_array.cc
                 blas_intfc23.cc
                 reorder_qt.cc
                 ci_wfn.cc
                 schmidt.cc
                 invert.cc
                 solve_pep.cc
                 #v_3.cc
                 david.cc
)
add_definitions("-DFC_SYMBOL=${FC_SYMBOL}")
psi4_add_module(lib qt sources_list psio ciomr mints)
