/*! \file
    \ingroup LMP2
    \brief Enter brief description of file here 
*/
#ifndef _psi_src_bin_lmp2_params_h_
#define _psi_src_bin_lmp2_params_h_
struct paramsinfo {
  int ref;
  char *wfn;
  char *jobtype;
  int print;
  int semicanonical;
  int opdm;
  int opdm_write;
  int opdm_print;
  int relax_opdm;
  long int memory;
  int nprocs;
  int myid;
  int ga_nprocs;
  int ga_myid;
  int num_threads;
  double rmsconv;
  double econv;
  int maxiter;
  int diis;
  int ndiis;
  int it_diis;
  double cutoff;
  double fskip;
};

#endif /* Header guard */
