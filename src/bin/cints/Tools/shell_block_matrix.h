#ifndef _psi_src_bin_cints_Tools_shell_block_matrix_h
#define _psi_src_bin_cints_Tools_shell_block_matrix_h

/*! \file shell_block_matrix.h
    \ingroup CINTS
*/

namespace psi { namespace CINTS {
double ****init_shell_block_matrix(void);
void free_shell_block_matrix(double****);
void shell_block_to_block(double****, double**);
void GplusGt(double****, double****);
};}
#endif
