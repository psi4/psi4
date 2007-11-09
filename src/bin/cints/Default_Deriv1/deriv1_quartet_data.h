#ifndef _psi_src_bin_cints_Default_Deriv1_deriv1_quartet_data_h
#define _psi_src_bin_cints_Default_Deriv1_deriv1_quartet_data_h

/*! \file deriv1_quartet_data.h
    \ingroup (CINTS)
    \brief Enter brief description of file here 
*/namespace psi { namespace CINTS {

void deriv1_quartet_data(prim_data *Data, double_array_t *fjt_table, double AB2, double CD2,
			 struct shell_pair *sp1, struct shell_pair *sp2, 
			 int am, int pi, int pj, int pk, int pl, double scale);
			 }}
			 #endif
