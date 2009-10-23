#ifndef _psi_src_bin_cints_MP2R12_r12_quartet_data_h
#define _psi_src_bin_cints_MP2R12_r12_quartet_data_h

/*! \file
    \ingroup CINTS
*/namespace psi { namespace cints {

void r12_quartet_data(prim_data* Data, double_array_t *fjt_table, double AB2, double CD2,
		      struct shell_pair* sp1, struct shell_pair* sp2, 
		      int am, int pi, int pj, int pk, int pl, double scale);
		      }}
		      #endif
