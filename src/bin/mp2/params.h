/*! \file
    \ingroup MP2
    \brief Enter brief description of file here 
*/
#ifndef _psi_src_bin_mp2_params_h_
#define _psi_src_bin_mp2_params_h_
struct params {
  int ref;
  char *wfn;
  char *jobtype;
  char *dertype;
  int cachelev;
  int cachetype;
  int print;
  int opdm;
  int opdm_write;
  int opdm_print;
  int relax_opdm;
  int gradient;
  long int memory;
  int semicanonical;
  int scs;
  double scs_scale_os;
  double scs_scale_ss;
};

#endif /* Header guard */
