/*! \file
    \ingroup CCENERGY
    \brief Enter brief description of file here 
*/

#ifndef _psi_src_bin_ccenergy_params_h
#define _psi_src_bin_ccenergy_params_h

namespace psi { namespace ccenergy {

/* Input parameters */
struct Params {
  int maxiter;
  double convergence;
  int restart;
  long int memory;
  char *aobasis;
  int cachelev;
  int cachetype;
  int ref;
  int diis;
  char *wfn;
  int print;
  int local;
  int num_amps;
  int print_mp2_amps;
  int brueckner;
  double bconv;
  int analyze;
  int print_pair_energies;
  int spinadapt_energies;
  int semicanonical;
  int local_mos;
  int dertype;
  int t2_coupled;
  char *prop;            /* user-selected property */
  int just_energy; /* just compute energy from T amplitudes on disk and quit */
	int just_residuals; /* just compute residuals from T amplitudes on disk and quit */
  char *abcd;
  int t3_Ws_incore;
  int nthreads;
  int scs;
  double scsmp2_scale_os;
  double scsmp2_scale_ss;
  double scscc_scale_os;
  double scscc_scale_ss;
  int newtrips;
};

}} // namespace psi::ccenergy

#endif // _psi_src_bin_ccenergy_params_h
