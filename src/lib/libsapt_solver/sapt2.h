#ifndef SAPT2_H
#define SAPT2_H

#include "sapt.h"

namespace boost {

template<class T> class shared_ptr;

}

namespace psi { namespace sapt {

class SAPT2 : public SAPT {
private:
  virtual void print_header();
  virtual void print_results();

  void df_integrals();
  void w_integrals();

protected:
  boost::shared_ptr<SAPTLaplaceDenominator> denom_;

  int nvec_;

  double **dAR_;
  double **dBS_;

  int *ioff_;
  int *index2i_;
  int *index2j_;

  int maxiter_;
  double e_conv_;
  double d_conv_;

  double e_elst10_;
  double e_exch10_;
  double e_exch10_s2_;
  double e_ind20_;
  double e_exch_ind20_;
  double e_disp20_;
  double e_exch_disp20_;
  double e_sapt0_;

  double **wBAA_;
  double **wBAR_;
  double **wBRR_;

  double **wABB_;
  double **wABS_;
  double **wASS_;

  double** get_AA_ints(const int);
  double** get_diag_AA_ints(const int);
  double** get_AR_ints(const int);
  double** get_RR_ints(const int);
  double** get_BB_ints(const int);
  double** get_diag_BB_ints(const int);
  double** get_BS_ints(const int);
  double** get_SS_ints(const int);
  double** get_AB_ints(const int);
  double** get_AS_ints(const int);
  double** get_RB_ints(const int);

  double **get_DF_ints(int, const char *, int, int, int, int);

  void cphf_solver(double**, double **, double *, int, const char *, 
    const char *, const char *, int, int);

  void exch_ind20rA_B();
  void exch_ind20rB_A();

  void tOVOV(int, const char *, int, int, int, double *, int, const char *, 
    int, int, int, double *, int, const char *);

public:
  SAPT2(Options& options, boost::shared_ptr<PSIO> psio, 
    boost::shared_ptr<Chkpt> chkpt);
  virtual ~SAPT2();

  virtual double compute_energy();

  void amplitudes();

  void elst10();
  void exch10_s2();
  void exch10();
  void ind20r();
  void exch_ind20r();
  void disp20();
  void exch_disp20();
};

}}

#endif
