/*! \file
    \ingroup TRANSQT
    \brief Enter brief description of file here 
*/
/* Struct for input parameters */

#ifndef _psi3_bin_transqt_Params_h_
#define _psi3_bin_transqt_Params_h_

namespace psi { namespace transqt {

struct Params {
  int runmode;                  /* which mode we are running in, defined
                                   in globals.h                           */
  int do_all_tei;               /* for two-elec ints, include fzc/fzv?    */
  int do_h_bare;                /* for fwd tf's, write bare h ints?       */
  int do_h_fzc;                 /* for fwd tf's, write dressed h fzc oper?*/
  int backtr;                   /* do a back-transformation?              */
  int do_oei;                   /* do one-electron integrals?             */
  int print_lvl;                /* defines verbosity of output            */
  int print_mos;                /* print SCF Molecular Orbitals?          */
  int print_te_ints;            /* print two-electron integrals?          */
  int print_oe_ints;            /* print one-electron integrals?          */
  int print_sorted_oe_ints;     /* print sorted one el integrals?         */
  int print_sorted_te_ints;     /* print sorted two el integrals?         */
  int max_buckets;              /* max number of Yoshimine bins           */
  int first_tmp_file;           /* number of first scratch binary file    */
  int presort_file;             /* filenumber for pre-sorted two el ints  */
  int keep_presort;             /* keep presort file after used?          */
  int jfile;                    /* filenumber for half-transformed ints   */
  int keep_half_tf;             /* keep file of half-transformed ints?    */
  int mfile;                    /* transformed two el ERIs                */
  int aa_mfile;                 /* transformed AA two el ERIs             */
  int bb_mfile;                 /* transformed BB two el ERIs             */
  int ab_mfile;                 /* transformed AB two el ERIs             */
  int h_bare_file;              /* filenum for bare onelec MO ints        */
  int h_bare_a_file;            /* filenum for alpha bare onelec MO ints  */
  int h_bare_b_file;            /* filenum for beta bare onelec MO ints   */
  int h_fzc_file;               /* file for dressed (fzc) onelec MO ints  */
  int h_fzc_a_file;             /* file for dressed (fzc) onelec MO ints  */
  int h_fzc_b_file;             /* file for dressed (fzc) onelec MO ints  */
  int sorted_tei_file;          /* filenumber for sorted two el ints      */
  /*    int src_ints_iwl;               use input ints in iwl format? Always! */
  int src_S_file;               /* filenumber for input overlap ints      */
  int src_T_file;               /* filenumber for input kinetic en. ints  */
  int src_V_file;               /* filenumber for input nuc. attr. ints   */
  int src_tei_file;             /* filenumber for input TEIs              */
  int tei_type;                 /* type of TEI integrals to do            */
  int tei_trans_type;           /* type of transformation to do           */
  int delete_src_oei;           /* delete input file of one el ints?      */
  int delete_src_tei;           /* delete input file of two el ints?      */
  int opdm_in_file;             /* filenum for MO one-particle dens mat   */
  int opdm_out_file;            /* filenum for AO one-particle dens mat   */
  int tpdm_add_ref;             /* add the reference contrib to the 2pdm? */
  int lag_in_file;              /* filenumber for Lagrangian (input)      */
  int lag_out_file;             /* filenumber for Lagrangian (output)     */
  char *wfn;                    /* string describing wavefunction type    */
  char *ref;                    /* string for reference type              */
  char *dertype;                /* derivative lvl: none, first, second..  */
  double tolerance;             /* tolerance on keeping integrals         */
  long int maxcor;              /* maximum available core memory in bytes */
  long int maxcord;             /* max ava core memory in doubles         */
  int fzc;                      /* really freeze core? (1 or 0)           */
  int del_restr_docc;           /* delete the restricted docc orbs? 1/0   */
  int treat_cor_as_fzc;         /* consider COR as FZC?  For DETCAS       */
  int print_reorder;            /* print the reordering array?            */
  int reorder;                  /* use user-given MO reordering array?    */
  int *moorder;                 /* user-given MO reordering array         */ 
  int lagran_double;            /* multiply the MO lagran by 2            */
  int lagran_halve;             /* multiply the MO lagran by 1/2          */
  int ras_type;                 /* define ras I to include 0 or excluded 1
				   socc                                   */
  int check_C_orthonorm;        /* check orthonormality of C matrix?      */
  int qrhf;                     /* boolean for QRHF reference             */
  int ivo;                      /* boolean for test IVO's                 */
  int semicanonical;            /* boolean for Semicanonical orbitals     */
                                /* now semicanonical is used as an integer:
                                  ==1 for semicanonical, 2 for z-averaged */

  char *aobasis;                /* string for AO-Basis CC algorithms      */

  int pitzer;                   /* boolean to override all MO reordering  */
  int psimrcc;                  /* boolean for psimrcc settings           */
};

/* Note that the current version does a reordering of orbital indices.
 * This is tantamount to sorting the one-electron intetrals.  However,
 * it is not quite "sorting" the two-electron integrals, so they are
 * currently output only to 'mfile' after reordering, and not to 
 * 'sorted_tei_file'
 */

}} // end namespace psi::transqt
#endif // header guard
