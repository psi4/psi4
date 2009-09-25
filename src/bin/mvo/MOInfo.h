/*! 
** \file
** \ingroup MVO
** \brief MOInfo structure of checkpoint information
*/

#ifndef _psi_src_bin_mvo_MOInfo_h
#define _psi_src_bin_mvo_MOInfo_h

namespace psi { namespace mvo {

/* Struct for PSIF_CHKPT molecular orbital information */
struct MOInfo {
  int nmo;               /* number of molecular orbitals                   */
  int nso;               /* number of basis functions (symmetry orbitals)  */
  int nao;               /* number of AO basis functions (nbfao)           */
  int nirreps;           /* number of irreducible reps                     */
  int iopen;             /* 1=open shell, 0=closed shell                   */
  int fzc_op_size;       /* size of frozen core operator                   */

  int nfzc;              /* number of frozen core orbitals, total          */
  int nfzv;              /* number of frozen virtual orbitals, total       */
  int ndocc;             /* number of doubly occupied orbitals, total      */
  int nsocc;             /* number of singly occupied orbitals, total      */

  int *sopi;             /* number of basis functions (SOs) per irrep      */
  int *orbspi;           /* number of molecular orbitals per irrep         */
  int *clsdpi;           /* number of closed-shell orbitals per irrep      */
  int *openpi;           /* number of open-shell orbitals per irrep        */
  int *virtpi;           /* number of virtual orbitals per irrep           */
  int *sosym;            /* array giving irrep for each SO                 */
  int *orbsym;           /* array giving irrep for each MO                 */

  int *order;            /* reordering array                               */
  int *corr2pitz_nofzv;  /* correlated->Pitzer order, excluding fzv's      */
  int *corr2pitz;        /* same as above but includes fzv's               */
  int *fruocc;           /* num of frozen virts per irrep                  */
  int *frdocc;           /* num of frozen core per irrep                   */
  int *rstrdocc;         /* num of restricted occupied per irrep           */
  int *rstruocc;         /* num of restricted unoccupied per irrep         */

  int *active;           /* num of active orbitals per irrep               */
  int *first;            /* first orbital (pitzer address) per irrep       */
  int *last;             /* last orbital per irrep                         */
  int *first_so;         /* first basis function (SO) per irrep            */
  int *last_so;          /* last basis function (SO) per irrep             */
  int *fstact;           /* first active orbital per irrep                 */
  int *lstact;           /* last active orbital per irrep                  */

  char **labels;         /* labels of the irreps                           */

  double ***evects;      /* SCF eigenvector matrix for each irrep          */
  double *fzc_operator;  /* AO frozen core operator (lwr triangle)         */
  int **ras_opi;         /* orbs per ras space per irrep                   */
};

}} // end namespace psi::mvo

#endif // header guard
