#ifndef SAPT2p_H
#define SAPT2p_H

#include "sapt2.h"

namespace boost {

template<class T> class shared_ptr;

}

namespace psi { namespace sapt {

class SAPT2p : public SAPT2 {
private:
  virtual void print_header();
  virtual void print_results();

protected:
  double e_disp21_;
  double e_disp22sdq_;
  double e_disp22t_;
  double e_est_disp22t_;
  double e_sapt2p_;

  void gARARxtARBS(int, const char *, const char, int, const char *, 
    const char *, const char *, int, int, int, int, int, int, int, 
    const char *);

  double disp21_1(int, const char *, const char *, int, int, int, int);
  double disp21_2(int, const char *, const char *, int, int);

  double disp211();
  double disp220s(int, const char *, const char *, int, const char *, 
    const char *, int, int, int);
  double disp220d_1(int, const char *, const char *, int, const char *,
    int, int, int);
  double disp220d_2(int, const char *, const char *, int, const char *,
    int, int, int, int, int, int, double *, double *, const char);
  double disp220q_1(int, const char *, const char *, const char *, int, int);
  double disp220q_2(int, const char *, const char *, const char *, int, 
    const char *, int, int, int);
  double disp220q_3(int, const char *, const char *, const char, int, 
    const char *, int, int, int, int, int, int);
  double disp220q_4(int, const char *, const char *, const char, int, 
    const char *, int, int, int, int, int, int);

  double disp220t(int, const char *, const char *, const char *, int, 
    const char *, int, const char *, int, int, int, int, int, int, double *, 
    double *);

public:
  SAPT2p(Options& options, boost::shared_ptr<PSIO> psio, 
    boost::shared_ptr<Chkpt> chkpt);
  virtual ~SAPT2p();

  virtual double compute_energy();

  virtual void amplitudes();

  void disp21();
  void disp22sdq();
  void disp22t();

};

}}

#endif
