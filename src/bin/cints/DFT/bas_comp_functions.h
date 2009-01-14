#ifndef _psi_src_bin_cints_DFT_bas_comp_functions_h
#define _psi_src_bin_cints_DFT_bas_comp_functions_h

/*! \file bas_comp_functions.h
    \ingroup CINTS
    \brief Enter brief description of file here 
*/
namespace psi { namespace CINTS {
  double calc_exp_basis(int shell_num, double rr);
  double calc_radial_bas(int shell_num, double rr, double r);
}}
#endif
