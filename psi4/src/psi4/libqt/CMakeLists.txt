set(sources_list lapack_intfc.cc
                 dx_write.cc
                 dirprd_block.cc
                 pople.cc
                 schmidt_add.cc
                 mat_print.cc
                 ras_set.cc
                 fill_sym_matrix.cc
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
                 schmidt.cc
                 invert.cc
                 solve_pep.cc
                 david.cc
)
add_definitions("-DFC_SYMBOL=${FC_SYMBOL}")
psi4_add_module(lib qt sources_list psio ciomr mints)
