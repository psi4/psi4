/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2016 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

/*! \file
    \ingroup DETCI
    \brief Enter brief description of file here
*/
/*
** STRUCTS.H
**
** C. David Sherrill
** Center for Computational Quantum Chemistry
** University of Georgia
** 1996
*/

#ifndef _psi_src_bin_detci_structs_h
#define _psi_src_bin_detci_structs_h

#include <string>
// #include <libmints/dimension.h>
#include <libmints/mints.h>

namespace psi { namespace detci {

/*** INCLUDES ***/
// do I really need these?  22 Jan 2008 CDS
// #include <unistd.h>
// #include <sys/<ctime>>

/*** DEFINES ***/
/*
typedef unsigned long long int BIGINT;
*/
typedef unsigned long int BIGINT;
//typedef boost::shared_ptr<psi::detci::CIvect> SharedCIVector;

#define CI_BLK_MAX 5000
#define IOFF_MAX 50604

#define PARM_GUESS_VEC_UNIT        0
#define PARM_GUESS_VEC_H0_BLOCK    1
#define PARM_GUESS_VEC_DFILE       3
#define PARM_OPENTYPE_UNKNOWN     -1
#define PARM_OPENTYPE_NONE         0
#define PARM_OPENTYPE_HIGHSPIN     1
#define PARM_OPENTYPE_SINGLET      2
#define PARM_OUTFILE_MAX           132
#define METHOD_RSP                 0
#define METHOD_OLSEN               1
#define METHOD_MITRUSHENKOV        2
#define METHOD_DAVIDSON_LIU_SEM    3
#define METHOD_RSPTEST_OF_SEM      4
#define PRECON_LANCZOS             0
#define PRECON_DAVIDSON            1
#define PRECON_EVANGELISTI         2
#define PRECON_GEN_DAVIDSON        3
#define PRECON_H0BLOCK_INVERT      4
#define PRECON_H0BLOCK_ITER_INVERT 5
#define PRECON_H0BLOCK_COUPLING    6
#define MULT                       0
#define DIV                        1
#define UPDATE_DAVIDSON            1
#define UPDATE_OLSEN               2
#define CI_VEC                     0
#define SIGMA_VEC                  1
#define TRUE                       1
#define FALSE                      0
#define HD_EXACT                   0
#define HD_KAVE                    1
#define ORB_ENER                   2
#define EVANGELISTI                3
#define LEININGER                  4
#define Z_HD_KAVE                  5

struct stringwr {
   unsigned char *occs;
   int **ij;
   int **oij;
   unsigned int **ridx;
   signed char **sgn;
   int *cnt;
   };


struct level {
   int num_j;
   int *a;
   int *b;
   int **k;
   int **kbar;
   int *y;
   int *x;
   };

struct stringgraph {
   int offset;
   int num_strings;
   struct level *lvl;
   int ***ktmp;             /* ktmp[case][row][level] */
   };


/*
** OLSEN_GRAPH structure.  This maintains a graphical
** representation of alpha and/or beta strings according to the method
** of Roos and Olsen, which maintains a different subgraph for each
** possible combination of point-group irrep, #RAS I electrons, and
** #RAS III electrons.
*/
struct olsen_graph {
   int num_str;             /* total number of strings */
   int num_drc_orbs;        /* number of dropped core orbitals */
   int num_expl_cor_orbs;   /* number of explicit core orbitals */
   int drc_sym;             /* symmetry (irrep) for the dropped core */
   int num_el;              /* number of electrons (total) in graph */
   int num_el_expl;         /* number of electrons (explicit) in graph */
   int num_orb;             /* number of orbitals explicitly treated */
   int ras1_lvl;            /* orbital number where RAS I ends (less drc),
                               or the last level in RAS I */
   int ras1_min;            /* minimum number of electrons in RAS I (for
                               the _strings_), incl. core */
   int ras1_max;            /* max number of RAS I electrons (useful when
                               the RAS I level may extend beyond the last
                               occupied orbital), incl. core */
   int ras3_lvl ;           /* orbital num where RAS III begins (less drc) */
   int ras3_max ;           /* maximum number of electrons in RAS III */
   int ras4_lvl ;           /* orbital number where RAS IV begins (less drc) */
   int ras4_max ;           /* maximum number of electrons in RAS IV */
   int nirreps;             /* number of irreps */
   int subgr_per_irrep;     /* possible number of Olsen subgraphs per irrep */
   int max_str_per_irrep;   /* largest number of strings found in an irrep */
   int *str_per_irrep;      /* array containing num strings per irrep */
   int ***decode;           /* decode[ras1_holes][ras3_e][ras4_e] */
   int **encode;            /* encode[0,1,2][code] gives ras1 e- (excl drc) and
                               ras3 e- and ras4 e- */
   struct stringgraph **sg; /* sg[irrep][code] */
   int *orbsym;             /* array for orbital irreps (incl. drc) */
   int *list_offset;        /* absolute offset for each list */
   };


/*
** FASTGRAPH structure.  Should provide a more straightforward (if
** also more restrictive) version of struct stringgraph.  Although
** somewhat more memory intensive, should allow for fast OTF
** computation of string addresses.  The memory requirements are not
** really severe in any case (0-5 MB total for all subgraphs).
*/
struct fastgraph {
   int num_strings;  /* number of strings in this subgraph */
   int **data;       /* holds k, then x, then y */
   };


/*
** GRAPH_SET structure.  This maintains a graphical
** representation of alpha and/or beta strings according to the method
** of Roos and Olsen, which maintains a different subgraph for each
** possible combination of point-group irrep, #RAS I electrons, and
** #RAS III electrons.  New implementation of previous olsen_graph
** structure; this one uses the more simplified fastgraph structure
** instead of stringgraph.
*/
struct graph_set {
   int num_str;             /* total number of strings */
   int num_graphs;          /* total number of valid subgraphs */
   int num_drc_orbs;        /* number of dropped core orbitals */
   int num_expl_cor_orbs;   /* number of explicit core orbitals */
   int drc_sym;             /* symmetry (irrep) for the dropped core */
   int num_el;              /* number of electrons (total) in graph */
   int num_el_expl;         /* number of electrons (explicit) in graph */
   int num_orb;             /* number of orbitals explicitly treated */
   int ras1_lvl;            /* orbital number where RAS I ends (not including
                               dropped core orbitals in the numbering),
                               or the last level in RAS I */
   int ras1_min;            /* minimum number of electrons in RAS I (for
                               the _strings_), no longer incl. core
                               (either explicit core or dropped core) */
   int ras1_max;            /* max number of RAS I electrons (useful when
                               the RAS I level may extend beyond the last
                               occupied orbital), no longer incl. core
                               (either explicit core or dropped core) */
   int ras3_lvl;            /* orbital num where RAS III begins (less drc) */
   int ras3_max;            /* maximum number of electrons in RAS III */
   int ras4_lvl;            /* orbital number where RAS IV begins (less drc) */
   int ras4_max;            /* maximum number of electrons in RAS IV */
   int ras34_max;           /* max number of electrons in RAS III AND IV */
   int nirreps;             /* number of irreps */
   int num_codes;           /* possible number of subgraphs per irrep */
   int ***decode;           /* decode[ras1_holes][ras3_e][ras4_e] */
   int **encode;            /* encode[0,1,2][code] gives ras1 e- (excl drc) and
                               ras3 e- and ras4 e- */
   struct fastgraph **AllGraph;
                            /* Pointers to all subgraphs */
   struct fastgraph **Graph;/* Pointers to allowed subgraphs */
   int *graph_irrep;        /* irrep of each non-null graph */
   int *graph_code;         /* code for each non-null graph */
   int *graph_offset;       /* absolute offset for each list */
   int *orbsym;             /* array for orbital irreps (incl. drc) */
   unsigned char ***Occs;   /* Orbital occupancies for each string */
   };




struct H_zero_block {
   double **H0b;               /* H0 block */
   double **H0b_inv;           /* inverse of block (H0 - E) */
   double **H0b_diag;          /* Eigenvectors of H0 block */
   double *H0b_diag_transpose; /* tmp array for Transpose of Eigenvectors of H0 block */
   double *H0b_eigvals;        /* Eigenvalues of H0 block */
   double *H00;                /* diag elements of H0 block */
   int size;                   /* size of H0 block */
   int osize;                  /* original (dimensioned size); can be reduced
                                  for Ms=0 cases by H0block_pairup() */
   int guess_size;             /* size of initial H0 block guess which may
                                  differ from size */
   int oguess_size;            /* original guess size may change with
                                  spin_cpl_chk() and pairup() */
   int coupling_size;          /* size of H0 block coupling to secondary
                                  space */
   int ocoupling_size;         /* original size of coupling_size */
   double *c0b, *s0b;          /* gathered C and sigma vectors */
   double *c0bp, *s0bp;        /* INV(H0 - E) times c0b and s0b */

   int *alplist;               /* list (graph) containing alpha string */
   int *betlist;               /* list containing beta string */
   int *alpidx;                /* relative index of alpha string */
   int *betidx;                /* relative index of beta string */
   int *blknum;                /* block number for each member */
   int *pair;                  /* which H0block member is related by
                                  interchange of alpha and beta indices */
   double **tmp1;              /* tmp matrix to hold (H0 - E) */
   int nbuf;                   /* number of buffers in CIvect */
   int *buf_num;               /* number of H0block elements per buffer */
   int **buf_member;           /* H0block members for each buffer */
   double spin_cp_vals;       /* Values of dets which should be added
                                  to the h0block but were not due to size
                                  restrictions of h0block.size */
   double *tmp_array1;         /* temporary array 1 */
   double *tmp_array2;         /* temporary array 2 */
   };


/*
** CalcInfo: Data Structure for holding calculation information such
**   as nuclear repulsion energy, number of atoms, number of basis functions,
**   etc.
*/
struct calcinfo {
   int natom;            /* number of atoms */
   int nso;              /* number of symmetry orbitals */
   int nmo;              /* number of molecular orbitals */
   int nmotri;           /* num elements in lwr diag matrix nmo big */
   int nirreps;          /* number of irreducible representations in pt grp */
   Dimension docc;            /* number of doubly occupied orbitals per irrep */
   Dimension socc;            /* number of singly occupied orbitals per irrep */
   Dimension ci_orbs;         /* number of act orbitals per irrep */
   Dimension dropped_docc;    /* number of core orbitals per irrep constrained to
                       be doubly occupied and dropped from explicit
                       consideration in the CI computation;
                       sum of frozen_docc and rstr_docc - CDS 4/15 */
   Dimension frozen_docc;     /* number of frozen doubly occupied orbs per irrep;
                       these are not explicitly present in the CI
                       and these orbitals are not allowed to optimize
                       in an MCSCF - CDS 4/15 */
   Dimension rstr_docc;       /* number of restricted doubly occupied orbs per
                       irrep; these are not explicitly present in the
                       CI and these orbitals are allowed to optimize
                       in an MCSCF - CDS 4/15 */
   Dimension dropped_uocc;     /* number of unoccupied orbitals per irrep
                       constrained to be unoccupied and dropped from
                       explicit consideratin in the CI computation;
                       sum of frozen_uocc and rstr_uocc - CDS 4/15 */
   Dimension frozen_uocc;     /* number of frozen unoccupied orbitals per irrep;
                       these are not explicitly present in the CI
                       and these orbitals are not allowed to optimize
                       in an MCSCF */
   Dimension rstr_uocc;       /* number of restricted unoccpied orbitals per irrep;
                            these are not explicitly present in the CI
                            and these orbitals are allowed to optimize in
                            an MCSCF */
   int iopen;            /* flag for whether open shell or not */
   double enuc;          /* nuclear repulsion energy */
   double escf;          /* scf energy */
   double efzc;          /* energy of the frozen core orbitals */
   double edrc;          /* energy of the dropped core orbitals */
   double e0;            /* E0, zeroth order energy */
   double e0_drc;        /* two times the sum of the dropped core orbitals */
   double e1;            /* E1, first order energy */
   int num_alp;          /* number of alpha electrons */
   int num_bet;          /* number of beta electrons */
   int num_alp_expl;     /* number of alpha electrons explicitly treated */
   int num_bet_expl;     /* number of beta electrons explicitly treated */
   char **labels;        /* labels for irreps */
   int *orbsym;          /* irrep for each orbital */
   std::vector<int> reorder;         /* map Pitzer-ordered orbitals to our ordering */
   std::vector<int> order;           /* map our ordering back to Pitzer ordering */
   std::vector<int> act_reorder;     /* map Pitzer-ordered orbitals to our ordering for active only*/
   std::vector<int> act_order;       /* map our ordering back to Pitzer ordering for active only*/
   std::vector<double> scfeigval;    /* SCF eigenvalues */
   std::vector<double> scfeigvala;    /* For ZAPTn, alpha and beta eigenvalues different */
   std::vector<double> scfeigvalb;    /* in SOCC space */
   SharedMatrix so_onel_ints; /* Pitzer-order one-electron integrals */
   SharedVector onel_ints;    /* CI-order one-electron integrals */
   SharedVector tf_onel_ints; /* CI-order transformed (avg) one-electron integrals */
   SharedVector gmat;         /* CI-order onel ints in RAS g matrix form, not symmetry packed */
   SharedVector twoel_ints;   /* CI-order two-electron integrals */
   int num_fzc_orbs;     /* number of frozen core orbitals */
   int num_rsc_orbs;     /* number of restricted core orbitals */
   int num_drc_orbs;     /* number of dropped core orbitals
                            (frozen + restricted) */
   int num_fzv_orbs;     /* number of frozen/deleted virtual orbitals */
   int num_rsv_orbs;     /* number of restricted virtual orbitals */
   int num_drv_orbs;     /* number of dropped virtual orbitals
                            (frozen + restricted) */
   int num_rot_orbs;     /* number of rotatable orbitals
                            (nmo - nfzc - nfzv) */
   int npop;             /* number of populated orbitals, nso - total virtual */
   int num_expl_cor_orbs;/* number of explicit core orbitals, i.e.,
                            orbitals that are constrained to be doubly
                            occupied in the CI and but are nevertheless
                            explicitly included in the CI computation (not
                            currently used for anything but kept because
                            lower-level routines allow it; it is set to 0) */
   int num_alp_str;      /* number of alpha strings */
   int num_bet_str;      /* number of beta strings */
   int num_ci_orbs;      /* nmo - num orbs frozen */
   size_t num_ci_tri;    /* number of upper triangular ci orbs */
   size_t num_ci_tri2;   /* number of doubly upper triangular ci orbs (example: two-el integrals) */
   int ref_alp;          /* address of reference alpha string */
   int ref_bet;          /* address of reference beta string */
   int ref_alp_list;     /* string list containing reference alpha string */
   int ref_bet_list;     /* string list containing reference beta string */
   int ref_alp_rel;      /* relative index of reference alpha string */
   int ref_bet_rel;      /* relative index of reference beta string */
   int ref_sym;          /* symmetry (irrep) of reference determinant */
   int spab;             /* socc per alpha or beta, for singlet states */
   int **ras_opi;        /* num orbs per irr per ras space ras_opi[ras][irr] */
   int **ras_orbs[4];    /* ras_orbs[ras][irr][cnt] gives an orbital number */
   int sigma_initialized; /* has sigma_init been called yet? */

};


/*
** parameters structure: holds run-time parameters
*/
struct params {
   std::string dertype;/* derivative level: none, first, etc */
   std::string wfn;    /* wavefunction type: CI, DETCAS, etc. */
   std::string ref;    /* reference type (RHF, ROHF); ROHF with MULTP=1
                          is an open-shell singlet */
   int multp;        /* multiplicity (2S+1) */
   int ex_lvl;       /* excitation level */
   int val_ex_lvl;   /* valence excitation level, used for RAS's */
   int cc_val_ex_lvl;/* NOT analogous to val_ex_lvl ... this controls how
                        many holes allowed from RAS II for MRCC */
   int cc_a_val_ex_lvl; /* alpha part of cc_val_ex_lvl */
   int cc_b_val_ex_lvl; /* beta part of cc_val_ex_lvl */
   int maxiter;      /* maximum number of allowed iterations */
   int num_roots;    /* number of CI roots to find */
   int istop;        /* stop after setting up CI space */
   int print_lvl;    /* print verbosity level */
   int print_ciblks; /* print a summary of the CI blocks? */
   bool die_if_not_converged;  /* Do we die if we do not converge */
   bool diag_h_converged;  /* Boolean to set if the diagonalization converges */
   double convergence;  /* convergence on RMS of the CI update vector */
                     /* (i.e. the Davidson/Liu d vector) applied to ea root */
   double energy_convergence;  /* convergence on CI energy */
   int ras;          /* do a RAS calculation?  Set true if "RAS1" keyword */
   int fci;          /* do a FULL ci calc?  (affects sigma1-2 subroutines) */
   int fci_strings;  /* do a FULL ci calc?  (affects string storage) */
   double S;         /* the value of quantum number S */
   int Ms0;          /* 1 if Ms=0, 0 otherwise */
   int ref_sym;      /* irrep for CI vectors;  -1 = find automatically */
   int opentype;     /* none, highspin, or open-shell singlet; see #define */
   int a_ras1_lvl;   /* orbital number defining RAS I for alpha electrons:
                        any electrons at this level or lower are in RAS I */
   int a_ras1_min;   /* minimum number of alpha electrons in RAS I */
   int a_ras1_max;   /* maximum number of alpha electrons in RAS I */
   int b_ras1_lvl;   /* orbital number defining RAS I for beta electrons:
                        any electrons at this level or lower are in RAS I */
   int b_ras1_min;   /* minimum number of beta electrons in RAS I */
   int b_ras1_max;   /* maximum number of beta electrons in RAS I */
   int a_ras3_max;   /* maximum number of alpha electrons in RAS III */
   int b_ras3_max;   /* maximum number of beta electrons in RAS III */
   int cc_a_ras3_max;/* as above but for CC */
   int cc_b_ras3_max;/* as above but for CC */
   int a_ras4_max;   /* maximum number of alpha electrons in RAS IV */
   int b_ras4_max;   /* maximum number of beta electrons in RAS IV */
   int cc_a_ras4_max;/* as above but for CC */
   int cc_b_ras4_max;/* as above but for CC */
   int ras1_lvl;     /* Orbital number of the highest orbital in RAS I */
   int ras1_min;     /* currently min #e AT THE RAS I LEVEL (excluding drc) */
   int ras3_lvl;     /* orbital number defining RAS III overall */
   int ras4_lvl;     /* orbital number defining RAS IV overall */
                     /* make larger than num_ci_orbs if there is no RAS IV */
   int ras3_max;     /* maximum number of electrons in RAS III */
   int cc_ras3_max;  /* as above but for CC */
   int ras4_max;     /* maximum number of electrons in RAS IV */
   int cc_ras4_max;  /* as above but for CC */
   int ras34_max;    /* max number of electrons in RAS III AND IV */
   int cc_ras34_max; /* as above but for CC */
   int a_ras34_max;  /* max number of alp electrons in RAS III AND IV */
   int b_ras34_max;  /* max number of bet electrons in RAS III AND IV */
   int cc_a_ras34_max;/* as above but for CC */
   int cc_b_ras34_max;/* as above but for CC */
   int guess_vector; /* what kind of CI vector to start with; see #define */
   int h0blocksize;  /* size of H0 block in preconditioner. */
   int h0guess_size; /* size of H0 block for initial guess */
   int h0block_coupling_size; /* size of coupling block in preconditioner */
   int h0block_coupling; /* 1 if true; 0 otherwise */
   int hd_ave;       /* how to average H diag energies over spin coupling
                        sets */
   int hd_otf;       /* 1 if diag energies computed on the fly 0 otherwise */
   int nodfile;      /* 1 if no dfile used 0 otherwise works for nroots=1 */
   int nprint;       /* number of important determinants to print out */
   int cc_nprint;    /* number of most important CC amps per ex lvl to print */
   int mixed;        /* 1=allow mixed excitations, 0=don't */
   int mixed4;       /* 1=allow mixed excitations in RAS IV, 0=don't */
   int repl_otf;     /* generate single-replacements on-the-fly */
   int r4s;          /* restrict strings with e- in RAS IV: i.e. if an
                        electron is in RAS IV, then the holes in RAS I
                        must equal the particles in RAS III + RAS IV else
                        the string is discarded */
   int calc_ssq;     /* calculate the value of <S^2> or not */
   int icore;        /* core option:
                            0 = RAS subblock at a time
                            1 = Entire CI vector at a time
                            2 = Symmetry block at a time */
   int diag_method;  /* diagonalization method:
                            0 = RSP
                            1 = Olsen
                            2 = Mitrushenkov
                            3 = Davidson/Liu SEM method
                            4 = Stupid test of Davidson/Liu SEM method */
   int precon;             /* preconditioner for diagonalization method
                              0 = Lanczos
                              1 = Davidson
                              2 = Generalized Davidson or H0block
                              3 = H0block inverse */
   int update;             /* update vector in diag method
                              1 = Davidson
                              2 = Olsen */
   int mcscf;             /* 1(0) if computing MCSCF is TRUE(FALSE) */
   int mpn;                /* 1(0) if computing mpn series is TRUE(FALSE) */
   int zaptn;              /* 1(0) if computing zaptn series is TRUE(FALSE) */
   int save_mpn2;          /* 0 = save MPn energy, 1 = save MP(2n-1) energy,
                              2 = save MP(2n-2) energy */
   int mpn_schmidt;        /* 1(0) if a orthonormal vector space is employed
                              rather than storing the kth order wfn */
   int wigner;             /* 1(0) if wigner formulas used in Empn series */
   int diag_iters_taken;   /* Number of diagonalization iterations taken */
   int maxnvect;           /* maximum number of b vectors for SEM method */
   int nunits;             /* num of tmp files to use for CI vects and such */
   int collapse_size;      /* how many vectors to collapse to in SEM */
   int lse_collapse;       /* iterations between lst sqr ext */
   int lse_iter;           /* iterations between last lst sqr ext */
   int lse;                /* 1(0) if lst sqr ext is TRUE or FALSE */
   double lse_tolerance;   /* energy converged to tol to perform lse */
   int neg_only;           /* 1(0) if get -(+) values of diag elements */
   int hd_filenum;  /* first tmp file for H diagonal */
   int c_filenum;   /* first tmp file for CI coeffs */
   int s_filenum;   /* first tmp file for sigma coeffs */
   int d_filenum;   /* first tmp file for D correction vectors */
   int num_init_vecs;      /* number of initial vectors for Davidson method */
   int restart;            /* restart flag, 0 or 1 */
   int bendazzoli;         /* use bendazzoli algorithm */
   int opdm;               /* call the opdm subroutine? */
   int opdm_diag;          /* get ci natural orbitals? */
   int opdm_ave;           /* average the opdm over several states */
   int transdens;          /* compute transition densities? */
   int dipmom;             /* compute dipole moment or transition dip mom?  */
   int tpdm;               /* call the tpdm subroutine? */
   int root;               /* which root to optimize (write opdm/tpdm for) */
   double perturbation_parameter; /* z in H = H0 + z * H1 */
   int z_scale_H;          /* 1(0) if pert. scaling used */
   int have_special_conv;  /* have a special convergence value from the
                              command line or the DETCASMAN driver? */
   double special_conv;    /* special convergence value */
   int nthreads;           /* number of threads to use in sigma routines */
   int sf_restrict;        /* 1 if restrict CI space (CI blocks) to
                              do only determinants (or their
                              spin-complements) in RASCI versions of
                              Krylov's SF CI */
   int print_sigma_overlap;/* Print sigma overlap matrix?  Test for Arteum */
   int filter_guess;       /* 1 if we want to filter out some of our guess
                              vectors by checking the phase of a pair of determinants */
   int filter_guess_sign;  /* the desired phase between dets 1 and 2 */
   int filter_guess_Ia;    /* absolute alpha string addr for determinant 1 */
   int filter_guess_Ib;    /* absolute beta string addr for determinant 1 */
   int filter_guess_Ja;    /* absolute alpha string addr for determinant 2 */
   int filter_guess_Jb;    /* absolute beta string addr for determinant 2 */
   int filter_guess_Iaridx;/* relative alpha string addr for det 1 */
   int filter_guess_Ibridx;/* relative beta string addr for det 1 */
   int filter_guess_Jaridx;/* relative alpha string addr for det 2 */
   int filter_guess_Jbridx;/* relative beta string addr for det 2 */
   int filter_guess_Iac;   /* string list number for alpha of det 1 */
   int filter_guess_Ibc;   /* string list number for beta of det 1 */
   int filter_guess_Jac;   /* string list number for alpha of det 2 */
   int filter_guess_Jbc;   /* string list number for beta of det 2 */
   int filter_guess_H0_det1; /* H0block determinant number for det 1 */
   int filter_guess_H0_det2; /* H0block determinant number for det 2 */
   int filter_zero_det;       /* zero out any particular determinant? */
   int filter_zero_det_Ia;    /* absolute alpha string addr for zero det */
   int filter_zero_det_Ib;    /* absolute beta  string addr for zero det */
   int filter_zero_det_Iac;   /* string list number for alpha of zero det */
   int filter_zero_det_Ibc;   /* string list number for beta of zero det */
   int filter_zero_det_Iaridx;/* relative alpha string for zero det */
   int filter_zero_det_Ibridx;/* relative beta string for zero det */
   int follow_vec_num;     /* num components in user-specified vec to follow */
   std::vector<double> follow_vec_coef;/* array of coefficients for vec to follow */
   std::vector<int> follow_vec_Ia;     /* array of absolute alpha strings for vector */
   std::vector<int> follow_vec_Ib;     /* array of absolute beta  strings for vector */
   std::vector<int> follow_vec_Iac;    /* array of alpha string lists for vector */
   std::vector<int> follow_vec_Ibc;    /* array of beta  string lists for vector */
   std::vector<int> follow_vec_Iaridx; /* array of alpha relative idx for vector */
   std::vector<int> follow_vec_Ibridx; /* array of beta  relative idx for vector */
   std::vector<int> ex_allow;          /* Determine nonstandard excitation types, such
                              as CID, CIST, CIDTQ, etc. Array is length
                              ex_lvl and each element is 1 or 0.  1 means
                              that excitation level is allowed. */
   std::vector<int> average_states;    /* which states to average in a SA calc */
   std::vector<double> average_weights;/* the weights for each state in a SA calc */
   int average_num;        /* length of the above two arrays */
   int cc;                 /* do coupled-cluster */
   int cc_ex_lvl;          /* coupled-cluster excitation level */
   int cc_export;          /* export a CC vector? */
   int cc_import;          /* import a CC vector? */
   int cc_fix_external;    /* fix amplitudes involving RAS I or IV ?        */
   int cc_fix_external_min;/* num external indices before amp gets fixed    */
   int cc_mixed;           /* ignore block if num holes in RAS I and II is
                              > cc_ex_lvl and if any indices correspond to
                              RAS I or IV (i.e., include only all-active
                              higher excitations)                           */
   int cc_update_eps;      /* update T amps with orb eigvals or not?        */
   int diis;               /* do DIIS?                                      */
   int diis_start;         /* how many diis vectors built up before start   */
   int diis_freq;          /* how many iters to go before a diis step       */
   int diis_min_vecs;      /* how many vectors required before do diis?     */
   int diis_max_vecs;      /* how many vectors maximum to hold?             */
   int cc_variational;     /* variational energy expression?                */

};

/*
** CI Vector structure which keeps track of how many
** symmetry or RAS blocks in a CI vector, the string graph
** codes for each block, offsets for each block (to compute absolute
** indices), and a decoder which takes two block codes (alpha and
** beta) and determines the CI vector block number.
*/
struct ci_blks {
    BIGINT vectlen;            /* total number of elements in the CI vector */
    int num_blocks;            /* number of blocks in the CI vector */
    int Ia_code[CI_BLK_MAX];   /* gives the block's alpha string code */
    int Ib_code[CI_BLK_MAX];   /* gives the block's beta string code */
    int Ia_size[CI_BLK_MAX];   /* num of alp strings in the block */
    int Ib_size[CI_BLK_MAX];   /* num of bet strings in the block */
    BIGINT offset[CI_BLK_MAX];  /* offset for absolute numbering */
    int **decode;              /* gives the block number for a given pair
                                   of alpha and beta codes */
    int num_alp_codes;         /* number of alpha codes in decode matrix */
    int num_bet_codes;         /* number of beta codes in decode matrix */
    int *first_iablk;          /* first blocknum for a given Ia irrep */
    int *last_iablk;           /* last blocknum for a given Ia irrep */
    int subgr_per_irrep;       /* possible number of Olsen subgraphs per irrep */
    int nirreps;               /* number of molecular irreps */
    int Ms0;                   /* 1 if Ms=0, 0 otherwise */
};

/*
** Struct to keep track of variable required in a sigma calculation
*/
struct sigma_data {
    double *F;
    int **Jij[2];
    int **Joij[2];
    int **Jridx[2];
    int *Jcnt[2];
    signed char **Jsgn[2];
    unsigned char **Toccs;
    double **transp_tmp;
    double **cprime;
    double **sprime;
    double *V, *Sgn;
    int *L, *R;
    int max_dim;
};

/*
** Structs for MCSCF variables
*/
struct mcscf_params {
    double rms_grad_convergence; /* convergence on RMS of orbital grad           */
    double energy_convergence;   /* convergence on CI energy                     */
    int max_iter;                /* maximum number of casscf iterations          */
    std::string mcscf_type;      /* Is this df?                                  */
    std::string algorithm;       /* What convergence algorithm do we use?        */
    double max_rot;              /* Maximum value in the rotation matrix         */

    bool orbital_so;             /* Do we do second-order orbital orbital?       */
    double so_start_grad;        /* RMS of orbital grad threshold for one-step   */
    double so_start_e;           /* energy convergence threshold for one-step    */

    int diis_start;              /* Start DIIS at this iteration                 */
    int diis_freq;               /* Do DIIS every n steps                        */
    int diis_max_vecs;           /* Maximum number of DIIS vectors               */
};

}} // namespace psi::detci

#endif // header guard
