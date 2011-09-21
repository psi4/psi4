/*! \file
    \ingroup ccresponse
    \brief Enter brief description of file here 
*/
#include <string>
namespace psi { namespace ccresponse {

struct Params {
  int print;             /* Output level control */
  long int memory;       /* Memory available (in bytes) */
  int cachelev;          /* cacheing level for libdpd */
  int ref;               /* reference determinant (0=RHF, 1=ROHF, 2=UHF) */
  double *omega;         /* energy of applied field (a.u) for dynamic polarizabilities */
  int nomega;            /* number of field energies desired */
  int maxiter;           /* maximum number of iterations allowed to converge perturbed amp eqns. */
  double convergence;    /* convergence criterion for perturbed wfns */
  int restart;           /* boolean for allowing a restart from on-disk amps */
  int diis;              /* boolean for using DIIS extrapolation */
  std::string prop;            /* user-selected property */
  int local;             /* boolean for simluation of local correlation */
  int analyze;
  int dertype;
  std::string gauge;           /* choice of gauge for optical rotation */
  std::string wfn;
  std::string abcd;
  int num_amps;
  int sekino;  /* Sekino-Bartlett size-extensive model-III */
  int linear;  /* Bartlett size-extensive (?) linear model */
};

}} // namespace psi::ccresponse
