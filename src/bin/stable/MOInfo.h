/*! 
** \file
** \ingroup STABLE
** \brief MOInfo structure with info about MO's from checkpoint
*/

namespace psi { namespace stable {

struct MOInfo {
  int nirreps;        /* no. of irreducible representations */
  int nmo;            /* no. of molecular orbitals */
  int nso;            /* no. of symmetry orbitals */
  int *orbspi;        /* no. of MOs per irrep */
  int *clsdpi;        /* no. of closed-shells per irrep  */
  int *openpi;        /* no. of open-shells per irrep */
  int *uoccpi;        /* no. of unoccupied orbitals per irrep  */
  int *frdocc;        /* no. of frozen core orbitals per irrep */
  int *fruocc;        /* no. of frozen unoccupied orbitals per irrep */
  char **labels;      /* irrep labels */
  int *occpi;         /* no. of occupied orbs. (incl. open) per irrep */
  int *aoccpi;        /* no. of alpha occupied orbs. (incl. open) per irrep */
  int *boccpi;        /* no. of beta occupied orbs. (incl. open) per irrep */
  int *virtpi;        /* no. of virtual orbs. (incl. open) per irrep */
  int *avirtpi;       /* no. of alpha virtual orbs. (incl. open) per irrep */
  int *bvirtpi;       /* no. of beta virtual orbs. (incl. open) per irrep */
  int *occ_sym;       /* relative occupied index symmetry */
  int *aocc_sym;      /* relative alpha occupied index symmetry */
  int *bocc_sym;      /* relative beta occupied index symmetry */
  int *vir_sym;       /* relative virtual index symmetry */
  int *avir_sym;      /* relative alpha virtual index symmetry */
  int *bvir_sym;      /* relative beta virtual index symmetry */
  int *occ_off;       /* occupied orbital offsets within each irrep */
  int *aocc_off;      /* occupied alpha orbital offsets within each irrep */
  int *bocc_off;      /* occupied beta orbital offsets within each irrep */
  int *vir_off;       /* virtual orbital offsets within each irrep */
  int *avir_off;      /* virtual alpha orbital offsets within each irrep */
  int *bvir_off;      /* virtual beta orbital offsets within each irrep */
  int *qt_occ;        /* CC->QT active occupied reordering array */
  int *qt_aocc;       /* CC->QT alpha active occupied reordering array */
  int *qt_bocc;       /* CC->QT beta active occupied reordering array */
  int *qt_vir;        /* CC->QT active virtiual reordering array */
  int *qt_avir;       /* CC->QT alpha active virtiual reordering array */
  int *qt_bvir;       /* CC->QT beta active virtiual reordering array */

  int *pitzer2qt;     /* Pitzer -> QT translation array */
  int *qt2pitzer;     /* QT -> Pitzer translation array */
                                                                                
  int *pitzer2qt_a;   /* Pitzer -> QT translation array for alpha orbitals */
  int *qt2pitzer_a;   /* QT -> Pitzer translation array for alpha orbitals */
  int *pitzer2qt_b;   /* Pitzer -> QT translation array for beta orbitals */
  int *qt2pitzer_b;   /* QT -> Pitzer translation array for beta orbitals */

  int *rank;          /* actual dimension of A in each irrep */
  double **A_evals;   /* lowest few eigenvalues of Hessian in each irrep */ 
  double **A_triplet_evals;  /* lowest few triplet eigenvalues of RHF Hessian 
                                in each irrep */
};

}} // namespace psi::stable
