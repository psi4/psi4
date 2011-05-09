#ifndef RPA_H
#define RPA_H

#include "cc.h"

namespace boost {
template<class T> class shared_ptr;
}

namespace psi {

class PSIO;
class Chkpt;

namespace dfcc {

class dRPA : public CC {
private:
  void print_header();

  double *tIAJB_;
  double *t2IAJB_;
  double *xIAJB_;

  double **B_p_IA_;
  double **Th_p_IA_;

protected:
  boost::shared_ptr<DFCCDIIS> diis_;

  double df_compute_energy();
  double cd_compute_energy();

  void apply_denom();

  double df_energy();
  double df_store_error_vecs();

public:
  dRPA(Options& options, boost::shared_ptr<PSIO> psio, boost::shared_ptr<Chkpt> chkpt);
  virtual ~dRPA();

  virtual double compute_energy();

};

}}

#endif
