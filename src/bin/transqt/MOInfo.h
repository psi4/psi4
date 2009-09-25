/*! 
** \file
** \ingroup TRANSQT
** \brief MOInfo structure
*/

#ifndef _psi3_bin_transqt_MOInfo_h_
#define _psi3_bin_transqt_MOInfo_h_

namespace psi { namespace transqt {

/* Struct for chkpt file molecular orbital information */
struct MOInfo {
  int nmo;               /* number of molecular orbitals                   */
  int nso;               /* number of basis functions (symmetry orbitals)  */
  int nao;               /* number of AO basis functions (nbfao)           */
  int nirreps;           /* number of irreducible reps                     */
  int iopen;             /* 1=open shell, 0=closed shell                   */
  int noeints;           /* num one electron ints (AO basis) or nbstri     */
  int nteints;           /* number of two electron ints (AO basis)         */

  int nfzc;              /* number of frozen core orbitals, total          */
  int nfzv;              /* number of frozen virtual orbitals, total       */
  int ndocc;             /* number of doubly occupied orbitals, total      */
  int nsocc;             /* number of singly occupied orbitals, total      */

  int nshell;            /* number of shells                               */
  int *sloc;             /* starting index for each shell of AOs           */
  int *snuc;             /* nuclear index for each shell                   */
  int *stype;            /* angular momentum quantum number for each shell */

  int *sopi;             /* number of basis functions (SOs) per irrep      */
  int *orbspi;           /* number of molecular orbitals per irrep         */
  int *clsdpi;           /* number of closed-shell orbitals per irrep      */
  int *openpi;           /* number of open-shell orbitals per irrep        */
  int *virtpi;           /* number of virtual orbitals per irrep           */
  int *sosym;            /* array giving irrep for each SO                 */
  int *orbsym;           /* array giving irrep for each MO                 */

  int *order;            /* reordering array                               */
  int *order_alpha;      /* alpha-orbital reordering array                 */
  int *order_beta;       /* beta-orbital reordering array                  */
  int *corr2pitz_nofzv;  /* correlated->Pitzer order, excluding fzv's      */
  int *corr2pitz;        /* same as above but includes fzv's               */
  int *corr2pitz_nofzv_a;/* alpha correlated->Pitzer order, excluding fzv's*/
  int *corr2pitz_a;      /* same as above but includes fzv's               */
  int *corr2pitz_nofzv_b;/* beta correlated->Pitzer order, excluding fzv's */
  int *corr2pitz_b;      /* same as above but includes fzv's               */
  int *frdocc;           /* num of frozen core per irrep                   */
  int *fruocc;           /* num of frozen virts per irrep                  */
  int *rstrdocc;         /* num of restricted occupied per irrep           */
  int *rstruocc;         /* num of restricted unoccupied per irrep         */
  int *active;           /* num of active orbitals per irrep               */
  int *first;            /* first orbital (pitzer address) per irrep       */
  int *last;             /* last orbital per irrep                         */
  int *first_alpha;      /* first alpha orbital (pitzer) per irrep         */
  int *last_alpha;       /* last alpha orbital per irrep                   */
  int *first_beta;       /* first beta orbital (pitzer) per irrep          */
  int *last_beta;        /* last beta orbital per irrep                    */
  int *first_so;         /* first basis function (SO) per irrep            */
  int *last_so;          /* last basis function (SO) per irrep             */
  int *fstact;           /* first active orbital per irrep                 */
  int *lstact;           /* last active orbital per irrep                  */

  int backtr_nirreps;    /* number of irreps to be used in the backtrans   */
  int *backtr_mo_first;  /* first mol orb per irrep for backtransforms     */
  int *backtr_mo_last;   /* last  mol orb per irrep for backtransforms     */
  int *backtr_mo_fstact; /* first active mol orb per irrep for backtrans   */
  int *backtr_mo_lstact; /* last  active mol orb per irrep for backtrans   */
  int *backtr_mo_orbspi; /* num mol orbs per irrep for backtransforms      */
  int *backtr_mo_active; /* num active mol orbs per irrep for backtr       */
  int *backtr_ao_first;  /* first ao orb per irrep for backtransforms      */
  int *backtr_ao_last;   /* last  ao orb per irrep for backtransforms      */
  int *backtr_ao_orbspi; /* num ao orbs per irrep for backtransforms       */
  int *backtr_ao_orbsym; /* orbital symmetry array for orbitals in ao bas  */

  char **labels;         /* labels of the irreps                           */

  double enuc;           /* nuclear repulsion energy                       */
  double escf;           /* SCF energy                                     */
  double efzc;           /* Frozen core energy                             */
  double ***evects;      /* SCF eigenvector matrix for each irrep          */
  double ***evects_alpha;/* SCF alpha eigenvector matrix for each irrep    */
  double ***evects_beta; /* SCF beta eigenvector matrix for each irrep     */
  double **scf_vector;   /* Full SCF eigenvector matrix                    */
  double **scf_vector_alpha; /* Full SCF eigenvector matrix                */
  double **scf_vector_beta; /* Full SCF eigenvector matrix                 */
  /* double *evals;          SCF eigenvalue array                           */
  double *oe_ints;       /* one-electron AO integrals (lwr triangle)       */
  double *S;             /* AO overlap matrix (lwr triangle)               */
  double *fzc_density;   /* AO frozen core density matrix (lwr triangle)   */
  double *fzc_density_alpha; /* AO frozen core density matrix (lwr tri)    */
  double *fzc_density_beta;   /* AO frozen core density matrix (lwr tri)   */ 
  double *fzc_operator;  /* AO frozen core operator (lwr triangle)         */
  double *fzc_operator_alpha;  /* AO frozen core operator (lwr tri)        */
  double *fzc_operator_beta;  /* AO frozen core operator (lwr tri)         */
  /* double *te_ints; */
};

}} // end namespace psi::transqt
#endif // header guard 
