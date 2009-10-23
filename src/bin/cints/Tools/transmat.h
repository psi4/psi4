#ifndef _psi_src_bin_cints_Tools_transmat_h
#define _psi_src_bin_cints_Tools_transmat_h

/*! \file transmat.h
    \ingroup CINTS
*/namespace psi { namespace cints {

double*** build_transmat(int *sym_oper, int nirreps, int max_am);
}}
#endif
