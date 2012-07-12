/*! \file
    \ingroup CCDENSITY
    \brief Enter brief description of file here 
*/

namespace psi { namespace ccdensity {

struct MOInfo {
    int nirreps;        /* no. of irreducible representations */
    int nmo;            /* no. of molecular orbitals */
    int nso;            /* no. of symmetry orbitals */
    int nactive;        /* no. of active orbitals */
    int iopen;          /* 0=closed shell; >0=open shell */
    int *orbspi;        /* no. of MOs per irrep */
    int *clsdpi;        /* no. of closed-shells per irrep excl. frdocc */
    int *openpi;        /* no. of open-shells per irrep */
    int *uoccpi;        /* no. of unoccupied orbitals per irrep excl. fruocc */
    int *frdocc;        /* no. of frozen core orbitals per irrep */
    int *fruocc;        /* no. of frozen unoccupied orbitals per irrep */
    char **labels;      /* irrep labels */
    int nfzc;           /* total no. of frozen core orbitals */
    int nfzv;           /* total no. of frozen virtual orbitals */
    int nclsd;          /* total no. of closd shells excl. frdocc */
    int nopen;          /* total no. of open shells  */
    int nuocc;          /* total no. of unoccupied shells excl. fruocc */
    int *occ_sym;       /* active occupied index symmetry */
    int *aocc_sym;      /* alpha active occupied index symmetry */
    int *bocc_sym;      /* beta active occupied index symmetry */
    int *vir_sym;       /* active virtual index symmetry */
    int *avir_sym;      /* alpha active virtual index symmetry */
    int *bvir_sym;      /* beta active virtual index symmetry */
    int sym;            /* symmetry of converged CCSD state */
    int *occpi;         /* no. of active occ. orbs. (incl. open) per irrep */
    int *aoccpi;        /* no. of alpha active occ. orbs. (incl. open) per irrep */
    int *boccpi;        /* no. of beta active occ. orbs. (incl. open) per irrep */
    int *virtpi;        /* no. of active virt. orbs. (incl. open) per irrep */
    int *avirtpi;       /* no. of alpha active virt. orbs. (incl. open) per irrep */
    int *bvirtpi;       /* no. of beta active virt. orbs. (incl. open) per irrep */
    int *occ_off;       /* occupied orbital offsets within each irrep */
    int *aocc_off;      /* alpha occupied orbital offsets within each irrep */
    int *bocc_off;      /* beta occupied orbital offsets within each irrep */
    int *vir_off;       /* virtual orbital offsets within each irrep */
    int *avir_off;      /* alpha virtual orbital offsets within each irrep */
    int *bvir_off;      /* beta virtual orbital offsets within each irrep */
    int *cc_occ;        /* QT->CC active occupied reordering array */
    int *cc_aocc;       /* QT->CC alpha active occupied reordering array */
    int *cc_bocc;       /* QT->CC beta active occupied reordering array */
    int *cc_vir;        /* QT->CC active virtiual reordering array */
    int *cc_avir;       /* QT->CC alpha active virtiual reordering array */
    int *cc_bvir;       /* QT->CC beta active virtiual reordering array */
    int *qt_occ;        /* CC->QT active occupied reordering array */
    int *qt_aocc;       /* CC->QT alpha active occupied reordering array */
    int *qt_bocc;       /* CC->QT beta active occupied reordering array */
    int *qt_vir;        /* CC->QT active virtiual reordering array */
    int *qt_avir;       /* CC->QT alpha active virtiual reordering array */
    int *qt_bvir;       /* CC->QT beta active virtiual reordering array */
    double enuc;        /* Nuclear repulsion energy */
    double escf;        /* SCF energy from chkpt */
    double eref;        /* Reference energy */
    double ecc;         /* CC energy (CC2, CCSD, or CC3) from ccenergy */
    double et;          /* (T) energy from cctriples */
    double **opdm;      /* Onepdm in the full (fzc+clsd+socc+uocc) space */
    double **opdm_a;    /* Alpha Onepdm in the full (fzc+clsd+socc+uocc) space */
    double **opdm_b;    /* Beta Onepdm in the full (fzc+clsd+socc+uocc) space */
    double **I;         /* Lagrangian matrix in the full space */
    double **I_a;       /* Alpha Lagrangian matrix in the full space */
    double **I_b;       /* Beta Lagrangian matrix in the full space */
  double **ltd;         /* <0|O|n> Left transition density */
  double **ltd_a;       /* <0|O|n> Left transition alpha density */
  double **ltd_b;       /* <0|O|n> Left transition beta density */
  double **rtd;         /* <n|O|0> Right transition density */ 
  double **rtd_a;       /* <n|O|0> Right transition alpha density */
  double **rtd_b;       /* <n|O|0> Right transition beta density */
  int *pitzer2qt;    /* Pitzer to QT re-ordering array */
  int *qt2pitzer;    /* QT to Pitzer re-ordering array */
  double **scf_qt;      /* SCF orbitals (QT ordering of MOs) */
  double ***L;
  double ***nabla;
  double ***dip;
};

}} // namespace psi::ccdensity
