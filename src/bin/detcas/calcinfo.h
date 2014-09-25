/*! \defgroup DETCAS detcas: Orbital optimizer for detci */

/*! 
  \file
  \ingroup DETCAS
  \brief Enter brief description of file here 
 
   CALCINFO.H
   
   C. David Sherrill
   University of California, Berkeley
   1998
*/

#ifndef _psi_src_bin_detcas_calcinfo_h
#define _psi_src_bin_detcas_calcinfo_h

namespace psi { namespace detcas {

/*
** CalcInfo: Data Structure for holding calculation information such
**   as nuclear repulsion energy, number of atoms, number of basis functions,
**   etc.
*/
struct calcinfo {
  int iter;              /* iteration number */
  int nmo;               /* number of molecular orbitals... the code often
                            uses nbfso instead but it shouldn't in case
                            of linear dependencies */
  int nbfso;             /* number of basis functions in symmetry orbitals */
  int nbstri;            /* num elements in lwr diag matrix nbfso big */
  int nbfao;             /* number of basis functions in atomic orbitals */
  int nirreps;           /* number of irreducible representations in pt grp */
  int num_fzc_orbs;      /* number of FZC orbitals (i.e. frozen core) */
  int num_cor_orbs;      /* number of COR orbitals (i.e. restricted core) */
  int num_vir_orbs;      /* number of VIR orbitals (i.e. restricted virtual) */
  int num_fzv_orbs;      /* number of frozen/deleted virtual orbitals */
  int npop;              /* number of populated orbitals, nbfso - nfzv */
  int max_orbs_per_irrep;/* max orbitals per irrep fzv not included */
  int max_pop_per_irrep; /* max populated orbitals per irrep fzv included */

  int *orbs_per_irr;     /* number of orbitals per irrep */
  int *docc;             /* doubly occupied orbitals per irrep */
  int *socc;             /* singly occupied orbitals per irrep */
  int *frozen_docc;      /* frozen doubly occupied orbs per irrep */
  int *frozen_uocc;      /* frozen virtual orbs per irrep */
  int *rstr_docc;        /* restricted doubly occupied orbs per irrep */
  int *rstr_uocc;        /* restricted virtual orbs per irrep */
  int *orbsym;           /* irrep for each orbital */
  int *pitz2ci;          /* map Pitzer-ordered orbitals to our ordering */
  int *ci2pitz;          /* map our ordering back to Pitzer ordering */
  int *ci2relpitz;       /* map CI ordering to _relative_ pitzer ordering */
  char **labels;         /* labels for irreps */
  int **ras_opi;         /* num orbs per irr per ras space ras_opi[ras][irr] */
  int **fzc_orbs;        /* frozen core orbitals numbers [irrep][orbnum] */
  int **cor_orbs;        /* restricted core orbitals numbers [irrep][orbnum] */
  int **vir_orbs;        /* restr virtual orbitals numbers [irrep][orbnum] */
  int **fzv_orbs;        /* frozen virtual orbitals numbers [irrep][orbnum] */

  int ***ras_orbs;       /* ras_orbs[ras][irr][cnt] gives an orbital number */

  int *first;            /* first orbital per irrep (in Pitzer order)    */
  int *last;             /* last  orbital per irrep (in Pitzer order)    */
  int *fstact;           /* first active orb per irrep (in Pitzer order) */
  int *lstact;           /* last  active orb per irrep (in Pitzer order) */
  int *active;           /* num active orbs per irrep                    */

  double enuc;           /* nuclear repulsion energy */
  double efzc;           /* frozen-core energy */
  double ***mo_coeffs;   /* matrix of molecular orbitals in Pitzer order */
  double *onel_ints;     /* one-electron integrals */
  double *twoel_ints;    /* two-electron integrals */
  double **opdm;         /* one-particle density matrix */
  double *tpdm;          /* two-particle density matrix */
  double **lag;          /* the MO Lagrangian */
  double *F_act;         /* Active Fock Matrix */
  double *mo_grad;       /* the MO gradient, dimension number of ind pairs */
  double mo_grad_rms;    /* the RMS of the MO gradient */
  double scaled_mo_grad_rms; 
  double *mo_hess_diag;  /* the MO Hessian, diagonal elements only         */
  double **mo_hess;      /* full MO Hessian */
  double *theta_cur;     /* current orbital rotation angles */
  double *theta_step;    /* step in orbital rotation angles */
  };

}} // end namespace psi::detcas

#endif // header guard

