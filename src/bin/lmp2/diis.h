/*! \file
    \ingroup LMP2
    \brief Enter brief description of file here 
*/
#ifndef _psi_src_bin_lmp2_diis_h_
#define _psi_src_bin_lmp2_diis_h_

struct diisinfo {

  int div;
//  int maxer;  maxer is going to be replaced with ndiis from params.h
  int matsize;
  int dmat1;
  int dmat2;
  int omat;
  int nmat;
  double ****error;
  double ****T_ext;

};

#endif /* Header guard */
