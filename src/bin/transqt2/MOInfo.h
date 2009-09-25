/*! 
** \file
** \ingroup TRANSQT2
** \brief MOInfo structure with information about MO's from checkpoint
*/

namespace psi {
  namespace transqt2 {

struct MOInfo {
  int nirreps;           /* no. of irreducible representations */
  int nmo;               /* no. of molecular orbitals */
  int nso;               /* no. of symmetry orbitals */
  int nao;               /* no. of atomic orbitals */
  int *sopi;             /* no. of SOs per irrep */
  int *mopi;             /* no. of MOs per irrep */
  int *sosym;            /* SO symmetry array */
  int *mosym;            /* MO symmetry array */
  int *actpi;            /* no. of active MOs per irrep */
  int *actsym;           /* active MO symmetry array */
  int *clsdpi;           /* no. of closed-shells per irrep */
  int *openpi;           /* no. of open-shells per irrep */
  int *uoccpi;           /* no. of unoccupied orbitals per irrep */
  int *frdocc;           /* no. of frozen core orbitals per irrep */
  int *fruocc;           /* no. of frozen unoccupied orbitals per irrep */
  int *core;             /* no. of "core" orbitals per irrep (for fzc op) */
  char **labels;         /* irrep labels */
  int nfzc;              /* total no. of frozen core orbitals */
  int ncore;             /* total no. of "core" orbs = frdocc + rstr_docc */
  int nfzv;              /* total no. of frozen virtual orbitals */
  int nactive;           /* no. of active MOs */

  double enuc;           /* Nuclear repulsion energy */
  double efzc;           /* Frozen core energy */

  int *pitz2corr_one;      /* one-electron integral reordering array (RHF): Pitzer MO -> corr */
  int *pitz2corr_one_A;    /* one-electron integral reordering array (UHF): Pitzer MO -> corr (alpha) */
  int *pitz2corr_one_B;    /* one-electron integral reordering array (UHF): Pitzer MO -> corr (beta) */
  int *pitz2corr_two;      /* two-electron integral reordering array (RHF): Pitzer MO -> corr */
  int *pitz2corr_two_A;    /* two-electron integral reordering array (UHF): Pitzer MO -> corr (alpha) */
  int *pitz2corr_two_B;    /* two-electron integral reordering array (UHF): Pitzer MO -> corr (beta) */

  double ***C;             /* Irrep-blocked MO/SO transform matrix (RHF) */
  double ***C_a;           /* Irrep-blocked alpha MO/SO transform matrix (UHF) */
  double ***C_b;           /* Irrep-blocked beta MO/SO transform matrix (UHF) */
  double **C_full;         /* Full MO/SO transform matrix (RHF) */
  double **C_full_a;       /* Full alpha MO/SO transform matrix (UHF) */
  double **C_full_b;       /* Full beta MO/SO transform matrix (UHF) */

  int *C_offset;           /* Column offset for skipping core orbs in transformations */
};

  } // namespace transqt2
} // namespace psi
