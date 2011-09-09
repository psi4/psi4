/*! \file
    \ingroup CCENERGY
    \brief Enter brief description of file here 
*/

#ifndef _psi_src_bin_ccenergy_params_h
#define _psi_src_bin_ccenergy_params_h

#include <string>

namespace psi { namespace ccenergy {

/* Input parameters */
struct Params {
  int maxiter;
  double convergence;
  int restart;
  long int memory;
  std::string aobasis;
  int cachelev;
  int cachetype;
  int ref;
  int diis;
  std::string wfn;
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
  std::string prop;            /* user-selected property */
  int just_energy; /* just compute energy from T amplitudes on disk and quit */
	int just_residuals; /* just compute residuals from T amplitudes on disk and quit */
  std::string abcd;
  int t3_Ws_incore;
  int nthreads;
  int scs;
  int scsn;
  int scscc;
  double scsmp2_scale_os;
  double scsmp2_scale_ss;
  double scscc_scale_os;
  double scscc_scale_ss;
  int newtrips;
};

}} // namespace psi::ccenergy

#endif // _psi_src_bin_ccenergy_params_h
