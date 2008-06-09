#ifndef _psi_src_bin_cints_Tools_transform_h
#define _psi_src_bin_cints_Tools_transform_h

/*! \file transform.h
    \ingroup CINTS
*/

namespace psi { namespace CINTS {
void transform_i(double *data, double *puream_data, double **c2p, int am_i, int nj, int nk, int nl);
void transform_j(double *data, double *puream_data, double **c2p, int am_j, int ni, int nk, int nl);
void transform_k(double *data, double *puream_data, double **c2p, int am_k, int ni, int nj, int nl);
void transform_l(double *data, double *puream_data, double **c2p, int am_l, int ni, int nj, int nk);
double *transform_ijkl(double *data, double *puream_data, double **cc2pp_ij, double **cc2pp_kl,
		       int trans_ij, int am_i, int am_j, int am_k, int am_l);
double *transform_ijkl_sparse(double *data, double *puream_data,
		       mat_elem **pp2cc_sparse, mat_elem **cc2pp_sparse,
		       int *pp2cc_rowlength, int *cc2pp_rowlength, int trans_ij,
		       int am_i, int am_j, int am_k, int am_l);
};}
#endif
