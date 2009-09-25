/*! \file
    \ingroup CCSORT
    \brief Enter brief description of file here 
*/

namespace psi { namespace ccsort {

struct MOInfo {
  int nirreps;           /* no. of irreducible representations */
  int nmo;               /* no. of molecular orbitals */
  int nso;               /* no. of symmetry orbitals */
  int nao;               /* no. of atomic orbitals */
  int noeints;           /* no. unique one-electron integrals (ex. fruocc) */
  int iopen;             /* 0=closed shell; >0=open shell */
  int *sopi;             /* no. of SOs per irrep */
  int *orbspi;           /* no. of MOs per irrep */
  int *clsdpi;           /* no. of closed-shells per irrep ex. frdocc */
  int *openpi;           /* no. of open-shells per irrep */
  int *uoccpi;           /* no. of unoccupied orbitals per irrep ex. fruocc */
  int *frdocc;           /* no. of frozen core orbitals per irrep */
  int *fruocc;           /* no. of frozen unoccupied orbitals per irrep */
  char **labels;         /* irrep labels */
  int nfzc;              /* total no. of frozen core orbitals */
  int nfzv;              /* total no. of frozen virtual orbitals */
  int nactive;           /* total no. of active orbitals */
  int *orbsym;           /* QT-ordered orbital symmetry array */

  double **PX;        /* MO-basis x-linear momentum ints (Pitzer order) */
  double **PY;        /* MO-basis y-linear momentum ints (Pitzer order) */
  double **PZ;        /* MO-basis z-linear momentum ints (Pitzer order) */
  int *pitzer2qt;     /* Pitzer -> QT reordering array */
  int *qt2pitzer;     /* QT -> Pitzer reordering array */

  int *pitz2qt;          /* RHF/ROHF Pitzer->QT translation lookup */
  int *pitz2qt_A;        /* UHF Alpha Pitzer->QT translation lookup */
  int *pitz2qt_B;        /* UHF Beta Pitzer->QT translation lookup */
  int *qt2pitz;          /* RHF QT->Pitzer translation lookup */
  int *qt2pitz_A;        /* UHF Alpha QT->Pitzer translation */
  int *qt2pitz_B;        /* UHF Beta QT->Pitzer translation */

  int *occ;              /* boolean array for active occ. orbs. */
  int *aocc;              /* boolean array for active occ. orbs. */
  int *bocc;              /* boolean array for active occ. orbs. */
  int *vir;              /* boolean array for active virt. orbs. */
  int *avir;              /* boolean array for active virt. orbs. */
  int *bvir;              /* boolean array for active virt. orbs. */
  int *socc;             /* boolean array for active socc. orbs. */
  int *all_occ;          /* boolean array for occ. orbs. (in the full space) */
  int *all_aocc;          /* boolean array for occ. orbs. (in the full space) */
  int *all_bocc;          /* boolean array for occ. orbs. (in the full space) */
  int *all_vir;          /* boolean array for virt. orbs. (in the full space) */
  int *all_avir;          /* boolean array for virt. orbs. (in the full space) */
  int *all_bvir;          /* boolean array for virt. orbs. (in the full space) */
  int *all_socc;         /* boolean array for socc. orbs. (in the full space) */
  int *frozen;           /* boolean array for frz. orbs (in the full space) */

  int *cc_occ;           /* QT->CC active occupied reordering array */
  int *cc_aocc;           /* QT->CC active occupied reordering array */
  int *cc_bocc;           /* QT->CC active occupied reordering array */
  int *cc_vir;           /* QT->CC active virtiual reordering array */
  int *cc_avir;           /* QT->CC active virtiual reordering array */
  int *cc_bvir;           /* QT->CC active virtiual reordering array */
  int *cc_allocc;        /* QT->CC all occupied reordering array */
  int *cc_allaocc;        /* QT->CC all occupied reordering array */
  int *cc_allbocc;        /* QT->CC all occupied reordering array */
  int *cc_allvir;        /* QT->CC all virtual reordering array */
  int *cc_allavir;        /* QT->CC all virtual reordering array */
  int *cc_allbvir;        /* QT->CC all virtual reordering array */

  int *qt_occ;           /* CC->QT active occupied reordering array */
  int *qt_aocc;           /* CC->QT active occupied reordering array */
  int *qt_bocc;           /* CC->QT active occupied reordering array */
  int *qt_vir;           /* CC->QT active virtiual reordering array */
  int *qt_avir;           /* CC->QT active virtiual reordering array */
  int *qt_bvir;           /* CC->QT active virtiual reordering array */
  int *qt_allocc;        /* CC->QT all occupied reordering array */
  int *qt_allaocc;        /* CC->QT all occupied reordering array */
  int *qt_allbocc;        /* CC->QT all occupied reordering array */
  int *qt_allvir;        /* CC->QT all virtual reordering array */
  int *qt_allavir;        /* CC->QT all virtual reordering array */
  int *qt_allbvir;        /* CC->QT all virtual reordering array */

  int *occ_sym;          /* CC active occupied index symmetry */
  int *aocc_sym;          /* CC active occupied index symmetry */
  int *bocc_sym;          /* CC active occupied index symmetry */
  int *vir_sym;          /* CC active virtual index symmetry */
  int *avir_sym;          /* CC active virtual index symmetry */
  int *bvir_sym;          /* CC active virtual index symmetry */
  int *allocc_sym;       /* CC all occupied index symmetry */
  int *allaocc_sym;       /* CC all occupied index symmetry */
  int *allbocc_sym;       /* CC all occupied index symmetry */
  int *allvir_sym;       /* CC all virtual index symmetry */
  int *allavir_sym;       /* CC all virtual index symmetry */
  int *allbvir_sym;       /* CC all virtual index symmetry */

  int *occpi;            /* no. of occupied orbs. (incl. open) per irrep */
  int *aoccpi;            /* no. of occupied orbs. (incl. open) per irrep */
  int *boccpi;            /* no. of occupied orbs. (incl. open) per irrep */
  int *virtpi;           /* no. of virtual orbs. (incl. open) per irrep */
  int *avirtpi;           /* no. of virtual orbs. (incl. open) per irrep */
  int *bvirtpi;           /* no. of virtual orbs. (incl. open) per irrep */
  int *all_occpi;        /* no. of occ. orbs. (incl. open and fzc) per irrep */
  int *all_aoccpi;        /* no. of occ. orbs. (incl. open and fzc) per irrep */
  int *all_boccpi;        /* no. of occ. orbs. (incl. open and fzc) per irrep */
  int *all_virtpi;       /* no. of virt. orbs. (incl. open and fzc) per irrep */
  int *all_avirtpi;       /* no. of virt. orbs. (incl. open and fzc) per irrep */
  int *all_bvirtpi;       /* no. of virt. orbs. (incl. open and fzc) per irrep */

  int *occ_off;          /* active occ. orbital offsets within each irrep */
  int *aocc_off;          /* active occ. orbital offsets within each irrep */
  int *bocc_off;          /* active occ. orbital offsets within each irrep */
  int *vir_off;         /* active virt. orbital offsets within each irrep */
  int *avir_off;         /* active virt. orbital offsets within each irrep */
  int *bvir_off;         /* active virt. orbital offsets within each irrep */
  int *all_occ_off;      /* all occ. orbital offsets within each irrep */
  int *all_aocc_off;      /* all occ. orbital offsets within each irrep */
  int *all_bocc_off;      /* all occ. orbital offsets within each irrep */
  int *all_vir_off;     /* all virt. orbital offsets within each irrep */
  int *all_avir_off;     /* all virt. orbital offsets within each irrep */
  int *all_bvir_off;     /* all virt. orbital offsets within each irrep */

  double enuc;           /* Nuclear repulsion energy */
  double efzc;           /* Frozen core energy */
  double eref;           /* The reference energy (computed here) */

  double **scf;
  double **usotao;
  double **MUX;
  double **MUY;
  double **MUZ;
  double **LX;
  double **LY;
  double **LZ;
  int irrep_x;
  int irrep_y;
  int irrep_z;
  int irrep_Rx;
  int irrep_Ry;
  int irrep_Rz;
  double ***C;           /* Virtual orbital transformation matrix (for AO-basis B terms) */
};

}} //namespace psi::ccsort
