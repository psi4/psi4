/*! \defgroup TRANSQT transqt: Integral Transformation Program */

/*!
** \file
** \ingroup TRANSQT
** \brief The SO-to-MO integral transformation program
**
** Program to transform one- and two-electron integrals over
** symmetry-orbitals to integrals over molecular orbitals.
**
** The original version of TRANSQT was based on the Saunders/van Lenthe
** algorithm, which should have been very fast, but adding the symmetry
** restrictions made the code so complicated that we couldn't figure
** out how to formulate the innermost loops as fast vector algorithms
** without a lot of trouble.
**
** This transformation is a straight "four single-index transformations"
** algorithm, modified to use matrix multiplications for each step.  Each half
** of the transformation is done with two consecutive matrix multiplications
** for each value of indices pq (or kl in the second half).  This actually
** does more work then the Elbert algorithm because it doesn't take full
** advantage of the permutational symmetry of the integrals, but should
** run faster on superscalar or vector processors with a well-written
** matrix multiplication algorithm.  We currently use mmult() written by Ed
** Seidl for the DEC3100.  It would be worthwhile to re-write these routines
** at a later date with more loops unrolled.  In addition, it should be
** possible to insert any highly-optimized matrix-multiplication routine
** in order to obtain peak performance.  Such a routine could even be
** FORTRAN-based, since the main arrays used here are constructed using
** block_matrix().
**
** This version uses symmetry simplifications for D2h and its subgroups,
** and runs out-of-core within the main loops.  The input list of integrals
** from file33 is presorted via a Yoshimine sort into all rs for a given
** pq.  This list is then half-transformed.  The resulting integrals are
** then sorted to all pq for a given kl.  This list is then transformed
** completely.  The final integrals are currently written as all ij for a
** given kl to file92 as an iwl buffer.
** TDC 5/17/95
**
** In addition, this version excludes frozen core and virtual orbitals
** from the transformation.  While this should cut back significantly
** on the computational expense, it does mean that changes to the active
** orbitals requires a complete re-transformation.
** TDC 9/19/95
**
** A restricted transformation which produces integrals of the type (ov|ov)
** where o=occupied and v=virtual has been included.  This is useful
** primarily for MP2 energy codes.  A separate routine, transform_two_mp2()
** has been added to simplify the structure of transform_two() somewhat.
** TDC 9/26/95
**
** Code updated with more even more of the options available in the
** Saunders-van Lenthe version of TRANSQT.  Output files are no longer
** hardwired, the one-electron integrals are output with the
** iwl_wrtone() function, implicit frozen core is treated properly
** (i.e. "FZC" orbitals), and restricted core transformations (i.e. "COR"
** orbitals) are also possible.
** CDS 10/6/95
**
** The COR orbitals stuff somehow got broken; fixed that, added the -quiet
** command-line option, and got the RAS reordering array in a much nicer
** way.
** CDS 6/26/96
**
** Finished a backtransformation routine that runs in C1 symmetry
** by folding the AO/SO transform into C.
** CDS 1/30/98
**
** More work done to clean up and generalize the code for forwards and
** backwards transformations.  Removed (or at least commented out)
** the evil hardwiring of moinfo.first, moinfo.last, etc, that I did
** last time for backtransforms.  Now there are backtransform-specific
** entries in moinfo.  Additionally, the code is able to handle all
** conceivable combinations of frozen core/frozen virtuals, and it
** knows to ignore these keywords and transform all integrals in the
** case of gradients or MCSCF type procedures.  Currently we assume
** that the one and twopdm's have frozen core indices but no frozen
** virt indices, whereas the Lagrangian has all indices.
** CDS 3/18/98
**
** Yet more work to streamline the one-electron integral
** transformations.  The code in transform_one() had gotten messy with
** so many possible cases.  Now we always run over all orbitals
** (including frozen core and virtuals) when writing out the MO
** one-electron ints.  Frozen orbitals can be filtered out when
** reading with the new iwl rdone() function.  We reserve one file for
** bare one-electron ints, and another for the frozen core operator
** (one-electron ints are dressed with two-elec ints involving frozen
** core orbitals).  Made file numbering for exported files use the new
** psifiles.h macros...no more hardwiring the numbers.
** CDS 4/29/98
**
** Added transformation routines necessary for MP2-R12 in standard
** approximation A. I decided to use command line options to pass
** necessary information to TRANSQT. See transqt.1 for an
** explanation. Had to add two new transform_two_ functions, to
** transform regular integrals, like ERIs and integrals over r12,
** which have bra-ket symmetry, and pesky [r12,T1] integrals which do
** not only have bra-ket symmetry but are also anti-symmetric with
** respect to a permutation of the first two indices. The required
** transformation is a partial 4-index transformation of OGOG type (2
** indices refer to occupied orbitals and two are general
** indices). Had to add a function to properly write
** half-(yosh_write_arr_mp2r12a) and fully-transformed integrals
** (iwl_buf_wrt_arr_mp2r12a in LIBIWL) to IWL buffers. Remember that
** the occupied indices are in QTS order and the general indices are
** in Pitzer order. Sad but true.
** EFV 08/05/99
**
** Added forward and backward transformations for UHF cases.
** TDC, 03
**
** Daniel Crawford, David Sherrill, and Justin Fermann
** Center for Computational Quantum Chemistry
** University of Georgia */

/* Max length of ioff array */
#define IOFF_MAX 32641

#include <cmath>

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <libpsio/psio.h>
#include <libciomr/libciomr.h>
#include <libchkpt/chkpt.h>
#include <libmints/wavefunction.h>
#include <libqt/qt.h>
#include <libiwl/iwl.h>
#include <psifiles.h>
#include "MOInfo.h"
#include "Params.h"
#include "globals.h"
#include "psi4-dec.h"

#include <libmints/wavefunction.h>
#include <libmints/basisset.h>

namespace psi { namespace transqt {

int *ioff;
struct MOInfo moinfo;
struct Params params;

/* averaged==1 for three-subspace canonicalized via averaged Fock
   for OPT1, OPT2 or ZAPT,
   averaged==0 for regular semicanonicalization for RMP, etc. */
void semicanonical_fock(int averaged);
void transform_one(void);
void transform_two(void);
void cleanup(void);
void fzc_density(int nirreps, int *frdocc, double *Pc, double **C,
  int *first, int *first_so, int *last_so, int *ioff);
void ivo_density(int nirreps, int *frdocc, int *docc, int *socc, double *P,
  double **C, int *first, int *first_so, int *last_so, int *ioff);
void transform_two_mp2(void);
void transform_two_mp2r12a_t(void);
void transform_two_mp2r12a_gr(void);

void init_io(Options & options);
void title(void);
void init_ioff(void);
void get_parameters(Options & options);
void print_parameters(void);
void get_moinfo(Options & options);
void get_one_electron_integrals(void);
void exit_io(void);
void get_reorder_array(Options& options);
double *** construct_evects(const char *spin, int nirreps, int *active, int *sopi,
  int *orbspi, int *first_so, int *last_so, int *first, int *last,
  int *fstact, int *lstact, int printflag);
void destruct_evects(int nirreps, double ***evects);
int check_C(int nso, int nmo, double **Cmat, double *S);

PsiReturnType transqt(Options & options)
{
  params.print_lvl = 1;
  init_io(options);
  title();
  init_ioff();
  get_parameters(options);
  print_parameters();
  get_moinfo(options);
  get_reorder_array(options);

  get_one_electron_integrals();

  if (params.check_C_orthonorm)
    check_C(moinfo.nso, moinfo.nmo, moinfo.scf_vector, moinfo.S);

  switch (params.tei_trans_type) {
  case MAKE_GGGG:
      transform_two();
      break;
  case MAKE_OVOV:
      transform_two_mp2();
      break;
  case MAKE_OGOG:
      if (params.tei_type == ERI ||
      params.tei_type == R12)
      transform_two_mp2r12a_gr();
      else if (params.tei_type == R12T1)
      transform_two_mp2r12a_t();
      break;
  }

  if (params.do_oei)
      transform_one();

  cleanup();
  exit_io();
  return Success;
}

void init_io(Options & options)
{
  params.print_lvl = options.get_int("PRINT_LVL"); // default is 1

  std::string tmpstring  = options.get_str("MODE"); // default TO_MO
  if (tmpstring == "TO_MO")
    params.backtr = 0;
  else
    params.backtr = 1;

  params.psimrcc = options.get_bool("PSIMRCC"); // default false

  params.runmode = MODE_NORMAL; // the default

  if (options["MP2R12A"].has_changed()) {
    // case 0: //--- transform ERIs ---
    if (tmpstring == "MP2R12AERI")
       params.runmode = MODE_MP2R12AERI;
    // case 1: //--- transform ints of r12
    else if (tmpstring == "MP2R12AR12")
      params.runmode = MODE_MP2R12AR12;
    // case 2: //--- transform ints of [r12,T1]
    else if (tmpstring == "MP2R12AR12T1")
      params.runmode = MODE_MP2R12AR12T1;
  }

  if (params.print_lvl > 0) tstart();
}

void title(void)
{
  if (params.print_lvl) {
    fprintf(outfile, "\n");
    fprintf(outfile,"\t**************************************************\n");
    fprintf(outfile,"\t* TRANSQT:  Program to transform integrals from  *\n");
    fprintf(outfile,"\t*           the SO basis to the MO basis.        *\n");
    fprintf(outfile,"\t*                                                *\n");
    fprintf(outfile,"\t*            Daniel, David, & Justin             *\n");
    fprintf(outfile,"\t*                   Sept  1995                   *\n");
    fprintf(outfile,"\t**************************************************\n");
    fprintf(outfile, "\n");
    }
}


void init_ioff(void)
{
  int i;
  ioff = (int *) malloc(IOFF_MAX * sizeof(int));
  if(ioff == NULL) {
    fprintf(stderr, "(transqt): error malloc'ing ioff array\n");
    abort();
  }
  ioff[0] = 0;
  for(i=1; i < IOFF_MAX; i++) {
    ioff[i] = ioff[i-1] + i;
  }
}


void exit_io(void)
{
  if (params.print_lvl) tstop();
}


void get_parameters(Options & options)
{
  int errcod;
  int tol, i;

  params.wfn = options.get_str("WFN");
  params.backtr = options.get_str("MODE") == "MO_TO_SO";
  params.dertype = options.get_str("DERTYPE");

  /* The defaults below depend on what run mode we are in */
  if (params.runmode == MODE_NORMAL) {
    /*The restricted transform does not currently work for the new
      MP2 module MLA June 27, 2003*/
    /*if (strcmp(params.wfn, "MP2") == 0 && strcmp(params.dertype, "NONE") == 0)
      params.tei_trans_type = MAKE_OVOV;
    else*/
    params.tei_trans_type = MAKE_GGGG;
    params.tei_type = ERI;
    params.src_tei_file = PSIF_SO_TEI;
    params.mfile = PSIF_MO_TEI;
    params.do_oei = 1;
  }
  else if (params.runmode == MODE_MP2R12AERI) {
    params.tei_trans_type = MAKE_OGOG;
    params.tei_type = ERI;
    params.src_tei_file = PSIF_SO_TEI;
    params.mfile = PSIF_MO_TEI;
    params.do_oei = 1;
  }
  else if (params.runmode == MODE_MP2R12AR12) {
    params.tei_trans_type = MAKE_OGOG;
    params.tei_type = ERI;
    params.src_tei_file = PSIF_SO_R12;
    params.mfile = PSIF_MO_R12;
    params.do_oei = 0;
  }
  else if (params.runmode == MODE_MP2R12AR12T1) {
    params.tei_trans_type = MAKE_OGOG;
    params.tei_type = R12T1;
    params.src_tei_file = PSIF_SO_R12T1;
    params.mfile = PSIF_MO_R12T2;
    params.do_oei = 0;
  }
  else {
    fprintf(outfile, "transqt: Unrecognized runmode %d\n", params.runmode);
    abort();
  }

  /* Adding for UHF capabilities -- TDC, 06/14/01 */
  params.ref = options.get_str("REFERENCE");

  params.semicanonical = 0;
  /* Semicanonical orbitals for perturbation theory */
  if(params.ref == "ROHF" && (params.wfn == "CCSD_T" || params.wfn == "MP2" ||
    params.wfn == "CC3" || params.wfn == "EOM_CC3" || params.wfn == "CC2" || params.wfn == "EOM_CC2")) {
    params.ref = "UHF";
    params.semicanonical = 1;
  }

  /* Averaged semicanonical orbitals for ZAPTn */
  if(params.ref == "ROHF" && params.wfn == "ZAPTN")
    params.semicanonical = 2;

  tol = options.get_int("TOLERANCE");
  params.tolerance = 1.0*pow(10.0,(double) -tol);

  params.h_bare_file = options.get_int("OEI_FILE");
  /* UHF additions, TDC 6/01 */
  params.h_bare_a_file = options.get_int("OEI_A_FILE");
  params.h_bare_b_file = options.get_int("OEI_B_FILE");

  params.h_fzc_file = options.get_int("FZC_FILE");
  /* UHF additions, TDC 6/01 */
  params.h_fzc_a_file = options.get_int("FZC_A_FILE");
  params.h_fzc_b_file = options.get_int("FZC_B_FILE");

  /* The sorted_tei values don't actually seem to be used */
  params.sorted_tei_file = options.get_int("SORTED_TEI_FILE");

  if (params.backtr) {
    params.src_tei_file = options.get_int("TPDM_FILE");
  }
  else {
    params.src_S_file   = options.get_int("SO_S_FILE");
    params.src_T_file   = options.get_int("SO_T_FILE");
    params.src_V_file   = options.get_int("SO_V_FILE");
    params.src_tei_file = options.get_int("SO_TEI_FILE");
  }

  params.first_tmp_file = options.get_int("FIRST_TMP_FILE");
  params.opdm_in_file   = options.get_int("OPDM_IN_FILE");
  params.opdm_out_file  = options.get_int("OPDM_OUT_FILE");
  params.lag_in_file    = options.get_int("LAG_IN_FILE");

  params.presort_file   = options.get_int("PRESORT_FILE");
  params.keep_presort   = options.get_bool("KEEP_PRESORT");
  params.jfile          = options.get_int("J_FILE");
  params.keep_half_tf   = options.get_bool("KEEP_J");

  if (params.backtr) params.mfile = PSIF_AO_TPDM;
  if (options["M_FILE"].has_changed())
    params.mfile = options.get_int("M_FILE");
  /* UHF additions, TDC, 6/01 */
  params.aa_mfile    = options.get_int("AA_M_FILE");
  params.bb_mfile    = options.get_int("BB_M_FILE");
  params.ab_mfile    = options.get_int("AB_M_FILE");

  params.max_buckets = options.get_int("MAX_BUCKETS");

  if ((params.wfn == "OOCCD" || params.wfn == "DETCAS" ||
       params.wfn == "BCCD"  || params.wfn == "BCCD_T" ||
       params.wfn == "CASSCF"|| params.wfn == "RASSCF") &&
      !params.backtr)
    params.delete_src_tei = 0;
  else
    params.delete_src_tei = 1;

  /* If AO-basis chosen, keep the SO_TEI file */
  params.aobasis = options.get_str("AO_BASIS"); // default NONE
  if(params.aobasis == "DISK")
    params.delete_src_tei = 0;

  if (!params.backtr) {
    if (options["DELETE_AO"].has_changed())
      params.delete_src_tei = options.get_bool("DELETE_AO");
  }
  else {
    if (options["DELETE_TPDM"].has_changed())
      params.delete_src_tei = options.get_bool("DELETE_TPDM");
  }

  params.print_te_ints   = options.get_bool("PRINT_TE_INTEGRALS");
  params.print_oe_ints   = options.get_bool("PRINT_OE_INTEGRALS");
  params.print_sorted_oe_ints = options.get_bool("PRINT_SORTED_OE_INTS");
  params.print_sorted_te_ints = options.get_bool("PRINT_SORTED_TE_INTS");
  params.print_mos       = options.get_bool("PRINT_MOS");

  if (params.wfn == "CI" || params.wfn == "DETCAS" || params.wfn == "CASSCF" ||
   params.wfn == "RASSCF" || params.wfn == "DETCI" || params.wfn == "ZAPTN") {
    params.lagran_double = 1;
    params.lagran_halve = 0;
    params.ras_type = 1;
  }
  else {
    params.lagran_double = 0;
    params.lagran_halve = 0;
    params.ras_type = 0;
  }
  if (options["LAGRAN_DOUBLE"].has_changed())
    params.lagran_double = options.get_bool("LAGRAN_DOUBLE");
  if (options["LAGRAN_HALVE"].has_changed())
    params.lagran_halve = options.get_bool("LAGRAN_HALVE");

  if ((params.wfn == "OOCCD" || params.dertype == "FIRST" || params.wfn == "DETCAS" ||
   params.wfn == "CASSCF" || params.wfn == "RASSCF") && !params.backtr)
    params.do_all_tei = 1;
  else params.do_all_tei = 0;
  if (params.runmode == MODE_MP2R12AERI || params.runmode == MODE_MP2R12AR12 ||
      params.runmode == MODE_MP2R12AR12T1)
    params.do_all_tei = 1;

  if (options["DO_ALL_TEI"].has_changed())
    params.do_all_tei = options.get_bool("DO_ALL_TEI");

  params.do_h_bare = 1;  params.do_h_fzc = 1;

  if (params.backtr) {
    params.do_h_bare = 0;
    params.do_h_fzc = 0;
  }
  else {
    if (params.wfn == "OOCCD") {
      params.do_h_bare = 1;
      params.do_h_fzc = 0;
    }
    else if ((params.wfn =="CI" || params.wfn == "DETCI" || params.wfn == "ZAPTN") &&
     params.dertype == "NONE") {
      params.do_h_bare = 0;
      params.do_h_fzc = 1;
    }
  }

  params.tpdm_add_ref = options.get_bool("TPDM_ADD_REF"); // default is false

  if (params.wfn =="DETCAS" || params.wfn == "CASSCF" || params.wfn == "RASSCF")
    params.treat_cor_as_fzc = 1;
  else
    params.treat_cor_as_fzc = 0;


  params.fzc = options.get_bool("FREEZE_CORE");
  if (params.backtr) params.fzc = 0; /* can't freeze core for backtr */

  /* Mark Hoffmann's code can't handle deleted core orbitals because
     his code wants them to be in the middle of the orbital ordering
     stack.  That's not going to work easily with TRANSQT unless we
     simply transform all of them.
  */
  if (params.wfn == "GVVPT2" || params.wfn == "MCSCF")
    params.fzc=0;

  /* remove restricted docc from RAS 1 ? */

  params.del_restr_docc = options.get_bool("DELETE_RESTR_DOCC"); // default is true
  if (params.backtr) params.del_restr_docc = 0;

  params.print_reorder = options.get_bool("PRINT_REORDER"); // default is false

  //fndcor(&(params.maxcor),infile,outfile);
  params.maxcor  = Process::environment.get_memory();
  params.maxcord = params.maxcor/sizeof(double);

  params.pitzer = options.get_bool("PITZER"); // default is false

  /* pitzer = true for SCF second derivative calculation */
  if(params.wfn == "SCF" && params.dertype == "SECOND")
    params.pitzer = 1;
  if(params.wfn == "SCF_MVD" && params.dertype == "FIRST")
    params.pitzer = 1;

  if(params.psimrcc)
    params.pitzer = 1;

  params.reorder = options.get_bool("REORDER"); //default is false

  params.check_C_orthonorm = options.get_bool("CHECK_C_ORTHONORM"); //default is false

  params.qrhf = options.get_bool("QRHF"); //default is false

  params.ivo  = options.get_bool("IVO"); //default is false

  return;
}


void print_parameters(void)
{

  if (params.print_lvl) {
      fprintf(outfile,"\tInput Parameters:\n");
      fprintf(outfile,"\t-----------------\n");
      fprintf(outfile,"\tWavefunction           =  %s\n", params.wfn.c_str());
      if(params.semicanonical==1) {
      fprintf(outfile,"\tReference orbitals     =  ROHF changed to UHF for Semicanonical Orbitals\n");
      }
      else {
      fprintf(outfile,"\tReference orbitals     =  %s\n", params.ref.c_str());
      }
      fprintf(outfile,"\tBacktrans              =  %s\n",
                               params.backtr ? "Yes" : "No");
      fprintf(outfile,"\tPrint MOs              =  %s\n",
                                  (params.print_mos ? "Yes": "No"));
      fprintf(outfile,"\tFreeze Core            =  %s\n",
                                  (params.fzc ? "Yes" : "No"));
      fprintf(outfile,"\tDelete Restricted Docc =  %s\n",
                          (params.del_restr_docc ? "Yes" : "No"));
      fprintf(outfile,"\tDo All TEI             =  %s\n",
                                  (params.do_all_tei ? "Yes" : "No"));
      fprintf(outfile,"\tMemory (Mbytes)        =  %5.1f\n",params.maxcor/1e6);
      fprintf(outfile,"\tMax Buckets            =  %d\n",params.max_buckets);
      fprintf(outfile,"\tFirst Tmp File         =  %d\n",params.first_tmp_file);
      fprintf(outfile,"\tPresort File           =  %d\n", params.presort_file);
      fprintf(outfile,"\tSource TEI File        =  %d\n",
                                  params.src_tei_file);
      fprintf(outfile,"\tOpdm In File           =  %d\n",params.opdm_in_file);
      fprintf(outfile,"\tOpdm Out File          =  %d\n",params.opdm_out_file);
      fprintf(outfile,"\tLag In File            =  %d\n", params.lag_in_file);
      fprintf(outfile,"\tKeep Presort           =  %s\n",
                                  (params.keep_presort ? "Yes" : "No"));
      fprintf(outfile,"\tJ File                 =  %d\n", params.jfile);
      fprintf(outfile,"\tKeep J                 =  %s\n",
                                  (params.keep_half_tf ? "Yes" : "No"));
      fprintf(outfile,"\tM File                 =  %d\n", params.mfile);
      fprintf(outfile,"\tBare OEI file          =  %d\n",
                                  params.h_bare_file);
      fprintf(outfile,"\tFrozen Core OEI file   =  %d\n",
                                  params.h_fzc_file);
      fprintf(outfile,"\tSorted TEI file        =  %d\n",
                                  params.sorted_tei_file);
      fprintf(outfile,"\tDelete TEI source file =  %s\n",
                                  (params.delete_src_tei ? "Yes" : "No"));
      fprintf(outfile,"\tAdd TPDM Ref Part      =  %s\n",
                                  (params.tpdm_add_ref ? "Yes" : "No"));
      fprintf(outfile,"\tDo Bare OEI tranform   =  %s\n",
                                  (params.do_h_bare ? "Yes" : "No"));
      fprintf(outfile,"\tDo FZC  OEI tranform   =  %s\n",
                                  (params.do_h_fzc ? "Yes" : "No"));
      fprintf(outfile,"\tTolerance              =  %3.1e\n", params.tolerance);
      fprintf(outfile,"\tPrint Level            =  %d\n", params.print_lvl);
      fprintf(outfile,"\tPrint TE Ints          =  %s\n",
                                  (params.print_te_ints ? "Yes": "No"));
      fprintf(outfile,"\tPrint OE Ints          =  %s\n",
                                  (params.print_oe_ints ? "Yes" : "No"));
      fprintf(outfile,"\tPrint Sorted TE Ints   =  %s\n",
                                  (params.print_sorted_te_ints ? "Yes": "No"));
      fprintf(outfile,"\tPrint Sorted OE Ints   =  %s\n",
                                  (params.print_sorted_oe_ints ? "Yes": "No"));
      fprintf(outfile,"\tReorder MOs            =  %s\n",
                                  (params.reorder ? "Yes" : "No"));
      fprintf(outfile,"\tCheck C Orthonormality =  %s\n",
                                  (params.check_C_orthonorm ? "Yes" : "No"));
      fprintf(outfile,"\tQRHF orbitals          =  %s\n",
                                  (params.qrhf ? "Yes" : "No"));
      fprintf(outfile,"\tIVO orbitals           =  %s\n",
                                  (params.ivo ? "Yes" : "No"));
      if(params.semicanonical) {
        if(params.semicanonical==1)
          fprintf(outfile,"\tSemicanonical orbitals =  Yes\n");
        else /* semicanonical == 2, hopefully */
          fprintf(outfile,"\tSemicanonical orbitals =  Z-Averaged\n");
      }
      fprintf(outfile,"\tPitzer                 =  %s\n",
                                  (params.pitzer ? "Yes" : "No"));
    }

  return;
}

void get_moinfo(Options & options)
{
  int i,j,k,h,errcod,size,row,col,p,q,offset,first_offset,last_offset,warned;
  int *tmpi, nopen;
  double **tmpmat, **so2ao;
  int **ras_opi;

  chkpt_init(PSIO_OPEN_OLD);
  moinfo.nmo = chkpt_rd_nmo();
  moinfo.nso = chkpt_rd_nso();
  moinfo.nao = chkpt_rd_nao();
  moinfo.nirreps = chkpt_rd_nirreps();
  moinfo.iopen = chkpt_rd_iopen();
  moinfo.labels = chkpt_rd_irr_labs();
  moinfo.sopi   = chkpt_rd_sopi();
  moinfo.orbspi = chkpt_rd_orbspi();
  moinfo.clsdpi = chkpt_rd_clsdpi();
  moinfo.openpi = chkpt_rd_openpi();
  moinfo.enuc = chkpt_rd_enuc();
  moinfo.escf = chkpt_rd_escf();
  moinfo.nshell = Process::environment.reference_wavefunction()->basisset()->nshell();
//  moinfo.sloc = Process::environment.reference_wavefunction()->basisset()->;
//  moinfo.snuc = chkpt_rd_snuc();
//  moinfo.stype = chkpt_rd_stype();
  moinfo.rstrdocc = init_int_array(moinfo.nirreps);
  moinfo.rstruocc = init_int_array(moinfo.nirreps);

  /* Needed for MO  reordering */
  if(params.ref == "UHF" && params.semicanonical == 0) {
    moinfo.scf_vector_alpha = chkpt_rd_alpha_scf();
    moinfo.scf_vector_beta = chkpt_rd_beta_scf();
  }
  else if(params.ref == "UHF" && params.semicanonical == 1){
    moinfo.scf_vector = chkpt_rd_scf();
  }
  else {
    moinfo.scf_vector = chkpt_rd_scf();
  }

  /* reorder the MOs if the user has requested it */
  if (params.reorder) {

    if(params.ref == "UHF" && params.semicanonical == 0) {
      fprintf(stderr, "ERROR: MO reordering not allowed for UHF references.\n");
      abort();
    }
    params.moorder = init_int_array(moinfo.nmo);
    if (!options["MOORDER"].has_changed() || options["MOORDER"].size() != moinfo.nmo)
      fprintf(outfile,"\nTrouble parsing MOORDER.  Ignoring it.\n");
    else {
      int *tmpv = options.get_int_array("MOORDER");
      memcpy(params.moorder, tmpv, moinfo.nmo*sizeof(int));
      delete tmpv;

      /* print the MOORDER array */
      fprintf(outfile, "\nMOORDER array: \n");
      for (i=0; i<moinfo.nmo; i++)
        fprintf(outfile, "%3d", params.moorder[i]);
      fprintf(outfile, "\n");

      /* change numbering in MOORDER from relative to absolute Pitzer */
      for (i=0,k=0,h=-1; i<moinfo.nirreps; i++) {
        for (j=0; j<moinfo.orbspi[i]; j++) {
          params.moorder[k++] += h;
        }
        h += j;
      }

      /* swap the rows of the SCF coefficient matrix */
      tmpmat = init_matrix(moinfo.nso,moinfo.nmo);
      for (i=0; i<moinfo.nso; i++)
        for (j=0; j<moinfo.nmo; j++)
          tmpmat[i][j] = moinfo.scf_vector[i][j];
      for (i=0; i<moinfo.nso; i++)
        for (j=0; j<moinfo.nmo; j++)
          moinfo.scf_vector[i][j] = tmpmat[i][params.moorder[j]];

      free_matrix(tmpmat, moinfo.nso);
    }
  }

  moinfo.sosym = init_int_array(moinfo.nso);
  for (i=0,k=0; i<moinfo.nirreps; i++) {
    for (j=0; j<moinfo.sopi[i]; j++,k++) {
      moinfo.sosym[k] = i;
    }
  }

  moinfo.orbsym = init_int_array(moinfo.nmo);
  for (i=0,k=0; i<moinfo.nirreps; i++) {
    for (j=0; j<moinfo.orbspi[i]; j++,k++) {
      moinfo.orbsym[k] = i;
    }
  }

  moinfo.frdocc = Process::environment.reference_wavefunction()->frzcpi();
  moinfo.fruocc = Process::environment.reference_wavefunction()->frzvpi();

  if (params.wfn == "CI" || params.wfn == "DETCI" || params.wfn == "GVVPT2"
      || params.wfn == "MCSCF" || params.wfn == "OOCCD" || params.wfn == "ZAPTN"
      || params.wfn == "CASSCF" || params.wfn == "RASSCF" || params.wfn == "DETCAS") {

    ras_opi = init_int_matrix(MAX_RAS_SPACES,moinfo.nirreps);
    tmpi = init_int_array(moinfo.nmo);

    if (params.wfn == "GVVPT2" || params.wfn == "MCSCF")
      i=1;
    else i=0;

    if (!ras_set2(moinfo.nirreps, moinfo.nmo, params.fzc,
                 params.del_restr_docc, moinfo.orbspi,
                 moinfo.clsdpi, moinfo.openpi, moinfo.frdocc, moinfo.fruocc,
                 moinfo.rstrdocc, moinfo.rstruocc, ras_opi, tmpi,
                 params.ras_type, i, options)) {
      //fprintf(outfile, "Error in ras_set().  Aborting.\n");
      //abort();
      throw PSIEXCEPTION("ras_set2() returned error code");
    }

    free_int_matrix(ras_opi);
    free(tmpi);

    if (params.treat_cor_as_fzc) {
      for (i=0; i<moinfo.nirreps; i++) {
        moinfo.frdocc[i] += moinfo.rstrdocc[i];
        moinfo.rstrdocc[i] = 0;
        if (params.backtr) {
          moinfo.fruocc[i] += moinfo.rstruocc[i];
          moinfo.rstruocc[i] = 0;
        }
      }
    }

  }

  if (!params.fzc) {
    for (i=0; i<moinfo.nirreps; i++) {
      moinfo.frdocc[i] = 0;
    }
  }

  // (ACS/FAE  more hacking for MRCC code...)
  if(params.psimrcc){
    params.fzc = 0;
    for (i=0; i<moinfo.nirreps; i++) {
      moinfo.frdocc[i] = moinfo.fruocc[i] = 0;
    }
  }

  moinfo.nfzc = 0;
  moinfo.nfzv = 0;
  for (i=0; i<moinfo.nirreps; i++) {
    moinfo.nfzc += moinfo.frdocc[i];
    moinfo.nfzv += moinfo.fruocc[i];
  }

  if ( chkpt_rd_override_occ() == 0 ) {
    moinfo.ndocc = 0;

    if (options["DOCC"].has_changed() && options["DOCC"].size() == moinfo.nirreps) {
      tmpi = options.get_int_array("DOCC");

      for (i=0,warned=0; i<moinfo.nirreps; i++) {
        if (tmpi[i] != moinfo.clsdpi[i] && !warned && !params.qrhf) {
          fprintf(outfile, "\tWarning: DOCC doesn't match chkpt file\n");
          warned = 1;
        }
        moinfo.clsdpi[i] = tmpi[i];
        moinfo.ndocc += tmpi[i];
      }
      delete tmpi;
    }

    if (options["SOCC"].has_changed() && options["SOCC"].size() == moinfo.nirreps) {
      tmpi = options.get_int_array("SOCC");

      for (i=0,warned=0; i<moinfo.nirreps; i++) {
        if (tmpi[i] != moinfo.openpi[i] && !warned && !params.qrhf) {
          fprintf(outfile, "\tWarning: SOCC doesn't match chkpt file\n");
          warned = 1;
        }
        moinfo.openpi[i] = tmpi[i];
        moinfo.nsocc += tmpi[i];
      }
      delete tmpi;
    }
  }

  /* Dump the new occupations to chkpt file if QRHF reference requested
     TDC,3/20/00 */
  if(params.qrhf) {
    chkpt_wt_clsdpi(moinfo.clsdpi);
    chkpt_wt_openpi(moinfo.openpi);
    for(i=0,nopen=0; i < moinfo.nirreps; i++) nopen += moinfo.openpi[i];
    chkpt_wt_iopen(nopen);
  }

  moinfo.virtpi = init_int_array(moinfo.nirreps);
  for(i=0; i < moinfo.nirreps; i++) {
    moinfo.virtpi[i] = moinfo.orbspi[i]-moinfo.clsdpi[i]-moinfo.openpi[i];
  }

  if (params.print_lvl) {
    fprintf(outfile,"\n\tChkpt File Parameters:\n");
    fprintf(outfile,"\t------------------\n");
    fprintf(outfile,"\tNumber of irreps = %d\n",moinfo.nirreps);
    fprintf(outfile,"\tNumber of SOs    = %d\n",moinfo.nso);
    fprintf(outfile,"\tNumber of MOs    = %d\n\n",moinfo.nmo);
    fprintf(outfile,
        "\tLabel\t# SOs\t# MOs\t# FZDC\t# DOCC\t# SOCC\t# VIRT\t# FZVR\n");
    fprintf(outfile,
        "\t-----\t-----\t-----\t------\t------\t------\t------\t------\n");
    for(i=0; i < moinfo.nirreps; i++) {
      fprintf(outfile,
          "\t %s\t   %d\t   %d\t    %d\t    %d\t    %d\t    %d\t    %d\n",
          moinfo.labels[i],moinfo.sopi[i],moinfo.orbspi[i],moinfo.frdocc[i],
          moinfo.clsdpi[i],moinfo.openpi[i],moinfo.virtpi[i],
          moinfo.fruocc[i]);
    }
    fprintf(outfile,"\n\tNuclear Repulsion Energy    = %20.10f\n",
        moinfo.enuc);
    fprintf(outfile,  "\tTotal SCF Energy            = %20.10f\n",
        moinfo.escf);
  }

  /* Set up some other useful parameters */
  moinfo.noeints = moinfo.nso*(moinfo.nso+1)/2;
  moinfo.nteints = moinfo.noeints*(moinfo.noeints+1)/2;

  if(params.semicanonical) {
    semicanonical_fock(params.semicanonical - 1); /* averaged == 0 for semicanonical, ==1 for Z-averaged canonicalization */
    if(params.semicanonical==1) {
        moinfo.scf_vector_alpha = chkpt_rd_alpha_scf();
        moinfo.scf_vector_beta = chkpt_rd_beta_scf();
    } else {
        moinfo.scf_vector = chkpt_rd_scf();
    }
  }

  /*
    Construct first and last index arrays for SOs: this defines the first
    absolute orbital index and last absolute orbital
    index for each irrep.  When there are no orbitals for an irrep, the
    value is -1 for first[] and -2 for last[].  Note that there must be
    basis functions in the first irrep (i.e. totally symmetric) for this to
    work.
  */
  moinfo.first_so = init_int_array(moinfo.nirreps);
  moinfo.last_so = init_int_array(moinfo.nirreps);
  for(h=0; h < moinfo.nirreps; h++) {
    moinfo.first_so[h] = -1;
    moinfo.last_so[h] = -2;
  }
  first_offset = 0;
  last_offset = moinfo.sopi[0] - 1;
  moinfo.first_so[0] = first_offset;
  moinfo.last_so[0] = last_offset;
  for(h=1; h < moinfo.nirreps; h++) {
    first_offset += moinfo.sopi[h-1];
    last_offset += moinfo.sopi[h];
    if(moinfo.sopi[h]) {
      moinfo.first_so[h] = first_offset;
      moinfo.last_so[h] = last_offset;
    }
  }

  /*
    Construct first and last index arrays: this defines the first
    absolute orbital index (Pitzer ordering) and last absolute orbital
    index for each irrep.  When there are no orbitals for an irrep, the
    value is -1 for first[] and -2 for last[].  Note that there must be
    orbitals in the first irrep (i.e. totally symmetric) for this to work.
  */
  moinfo.first = init_int_array(moinfo.nirreps);
  moinfo.last = init_int_array(moinfo.nirreps);
  for(h=0; h < moinfo.nirreps; h++) {
    moinfo.first[h] = -1;
    moinfo.last[h] = -2;
  }
  first_offset = 0;
  last_offset = moinfo.orbspi[0] - 1;
  moinfo.first[0] = first_offset;
  moinfo.last[0] = last_offset;
  for(h=1; h < moinfo.nirreps; h++) {
    first_offset += moinfo.orbspi[h-1];
    last_offset += moinfo.orbspi[h];
    if(moinfo.orbspi[h]) {
      moinfo.first[h] = first_offset;
      moinfo.last[h] = last_offset;
    }
  }
  /*
    Construct first and last active index arrays: this defines the first
    absolute orbital index (Pitzer ordering) and last absolute orbital
    index for each irrep, excluding frozen orbitals.  When there are no
    orbitals for an irrep, the value is -1 for first[] and -2 for last[].
    Note that there must be orbitals in the first irrep (i.e. totally
    symmetric) for this to work.
  */
  moinfo.fstact = init_int_array(moinfo.nirreps);
  moinfo.lstact = init_int_array(moinfo.nirreps);
  for(h=0; h < moinfo.nirreps; h++) {
    moinfo.fstact[h] = -1;
    moinfo.lstact[h] = -2;
  }
  first_offset = moinfo.frdocc[0];
  last_offset = moinfo.orbspi[0] - moinfo.fruocc[0] - 1;
  moinfo.fstact[0] = first_offset;
  moinfo.lstact[0] = last_offset;
  for(h=1; h < moinfo.nirreps; h++) {
    first_offset += moinfo.orbspi[h-1]+moinfo.frdocc[h]-moinfo.frdocc[h-1];
    last_offset += moinfo.orbspi[h] - moinfo.fruocc[h] + moinfo.fruocc[h-1];
    if(moinfo.orbspi[h]) {
      moinfo.fstact[h] = first_offset;
      moinfo.lstact[h] = last_offset;
    }
  }

  /* Now define active[] such that frozen orbitals are taken into account */
  moinfo.active = init_int_array(moinfo.nirreps);
  for(h=0; h < moinfo.nirreps; h++) {
    moinfo.active[h] = moinfo.orbspi[h]-moinfo.frdocc[h]-moinfo.fruocc[h];
  }

  /* Organize MOs into symmetry blocks; good for symmetry-blocked transform */
  /* We need to really change things around at this point if it is a        */
  /*   backtransformation  -- run everything in effectively C1 symmetry     */
  /*   and put the SO/AO transformation into the SCF coefficient matrix.    */
  /* It is possible to maintain symmetry by doing an MO->SO->AO backtrans,  */
  /*   but this requires some complicated sorting that we will avoid now.   */

  if(params.print_mos) {
    fprintf(outfile,"\n\tSCF Eigenvectors (Transformation Matrix):\n");
    fprintf(outfile,  "\t-----------------------------------------\n");
  }

  if (!params.backtr) {
    if (params.do_all_tei) {

      if(params.ref == "UHF") {
    moinfo.evects_alpha = construct_evects("alpha", moinfo.nirreps, moinfo.orbspi,
                           moinfo.sopi, moinfo.orbspi,
                           moinfo.first_so, moinfo.last_so,
                           moinfo.first, moinfo.last,
                           moinfo.first, moinfo.last, params.print_mos);

    moinfo.evects_beta = construct_evects("beta", moinfo.nirreps, moinfo.orbspi,
                          moinfo.sopi, moinfo.orbspi,
                          moinfo.first_so, moinfo.last_so,
                          moinfo.first, moinfo.last,
                          moinfo.first, moinfo.last, params.print_mos);
      }
      else {
    moinfo.evects = construct_evects("RHF", moinfo.nirreps, moinfo.orbspi,
                     moinfo.sopi, moinfo.orbspi,
                     moinfo.first_so, moinfo.last_so,
                     moinfo.first, moinfo.last,
                     moinfo.first, moinfo.last, params.print_mos);

      }
    }
    else {
      if(params.ref == "UHF") {
    moinfo.evects_alpha = construct_evects("alpha", moinfo.nirreps, moinfo.active,
                           moinfo.sopi, moinfo.orbspi,
                           moinfo.first_so, moinfo.last_so,
                           moinfo.first, moinfo.last,
                           moinfo.fstact, moinfo.lstact, params.print_mos);
    moinfo.evects_beta = construct_evects("beta", moinfo.nirreps, moinfo.active,
                          moinfo.sopi, moinfo.orbspi,
                          moinfo.first_so, moinfo.last_so,
                          moinfo.first, moinfo.last,
                          moinfo.fstact, moinfo.lstact, params.print_mos);
      }
      else {
    moinfo.evects = construct_evects("RHF", moinfo.nirreps, moinfo.active,
                     moinfo.sopi, moinfo.orbspi,
                     moinfo.first_so, moinfo.last_so,
                     moinfo.first, moinfo.last,
                     moinfo.fstact, moinfo.lstact, params.print_mos);

      }
    }

  }

  /* For a backtransform, we need to fill the SCF coefficient matrix with
   *  the MO->AO matrix....this matrix is NOT symmetry-blocked.  We assume
   *  that the integrals coming in have been renumbered such that frozen
   *  virtual orbitals are all up at the top and that these orbitals make
   *  no contributions (those labels will not appear on input OPDM or TPDM).
   */
  else {

    if (params.print_mos) fprintf(outfile, "SO to AO matrix\n");
    so2ao = chkpt_rd_usotao();
    if (params.print_mos) print_mat(so2ao,moinfo.nso,moinfo.nao,outfile);
    tmpmat = init_matrix(moinfo.nso, moinfo.nmo - moinfo.nfzv);

    if(params.ref == "UHF") {
      moinfo.evects_alpha = (double ***) malloc (1 * sizeof(double **));
      moinfo.evects_alpha[0] = block_matrix(moinfo.nao, moinfo.nmo - moinfo.nfzv);
      moinfo.evects_beta = (double ***) malloc (1 * sizeof(double **));
      moinfo.evects_beta[0] = block_matrix(moinfo.nao, moinfo.nmo - moinfo.nfzv);

      /*** alpha SCF matrix ***/

      /* fill up a temporary SCF matrix with frozen virt columns deleted */
      for (h=0,offset=0; h < moinfo.nirreps; h++) {
    if (h > 0) offset += moinfo.fruocc[h-1];
    if (moinfo.first[h] < 0 || moinfo.lstact[h] < 0) continue;
    for (p=moinfo.first_so[h]; p <= moinfo.last_so[h]; p++) {
      for (q=moinfo.first[h]; q <= moinfo.lstact[h]; q++) {
        tmpmat[p][q-offset] = moinfo.scf_vector_alpha[p][q];
      }
    }
      }

      /* now that we have C, multiply it by the SO->AO transform matrix */
      mmult(so2ao,1,tmpmat,0,moinfo.evects_alpha[0],0,moinfo.nao,moinfo.nso,
        moinfo.nmo - moinfo.nfzv,0);
      if (params.print_mos) {
    fprintf(outfile, "Alpha C matrix (including AO to SO)\n");
    print_mat(moinfo.evects_alpha[0],moinfo.nao,moinfo.nmo-moinfo.nfzv,outfile);
      }

      /*** beta SCF matrix ***/

      zero_mat(tmpmat, moinfo.nso, moinfo.nmo - moinfo.nfzv);

      /* fill up a temporary SCF matrix with frozen virt columns deleted */
      for (h=0,offset=0; h < moinfo.nirreps; h++) {
    if (h > 0) offset += moinfo.fruocc[h-1];
    if (moinfo.first[h] < 0 || moinfo.lstact[h] < 0) continue;
    for (p=moinfo.first_so[h]; p <= moinfo.last_so[h]; p++) {
      for (q=moinfo.first[h]; q <= moinfo.lstact[h]; q++) {
        tmpmat[p][q-offset] = moinfo.scf_vector_beta[p][q];
      }
    }
      }

      /* now that we have C, multiply it by the SO->AO transform matrix */
      mmult(so2ao,1,tmpmat,0,moinfo.evects_beta[0],0,moinfo.nao,moinfo.nso,
        moinfo.nmo - moinfo.nfzv,0);
      if (params.print_mos) {
    fprintf(outfile, "Beta C matrix (including AO to SO)\n");
    print_mat(moinfo.evects_beta[0],moinfo.nao,moinfo.nmo-moinfo.nfzv,outfile);
      }

    }
    else {
      moinfo.evects = (double ***) malloc (1 * sizeof(double **));
      moinfo.evects[0] = block_matrix(moinfo.nao, moinfo.nmo - moinfo.nfzv);

      /* fill up a temporary SCF matrix with frozen virt columns deleted */
      for (h=0,offset=0; h < moinfo.nirreps; h++) {
    if (h > 0) offset += moinfo.fruocc[h-1];
    if (moinfo.first[h] < 0 || moinfo.lstact[h] < 0) continue;
    for (p=moinfo.first_so[h]; p <= moinfo.last_so[h]; p++) {
      for (q=moinfo.first[h]; q <= moinfo.lstact[h]; q++) {
        tmpmat[p][q-offset] = moinfo.scf_vector[p][q];
      }
    }
      }

      /* now that we have C, multiply it by the SO->AO transform matrix */
      mmult(so2ao,1,tmpmat,0,moinfo.evects[0],0,moinfo.nao,moinfo.nso,
        moinfo.nmo - moinfo.nfzv,0);
      if (params.print_mos) {
    fprintf(outfile, "C matrix (including AO to SO)\n");
    print_mat(moinfo.evects[0],moinfo.nao,moinfo.nmo-moinfo.nfzv,outfile);
      }
    }

    free_matrix(tmpmat, moinfo.nso);
    free_block(so2ao);

  }

  chkpt_close();

  /* define some first/last arrays for backtransforms */
  moinfo.backtr_nirreps = 1;
  moinfo.backtr_mo_first = init_int_array(moinfo.backtr_nirreps);
  moinfo.backtr_mo_last = init_int_array(moinfo.backtr_nirreps);
  moinfo.backtr_mo_fstact = init_int_array(moinfo.backtr_nirreps);
  moinfo.backtr_mo_lstact = init_int_array(moinfo.backtr_nirreps);
  moinfo.backtr_mo_orbspi = init_int_array(moinfo.backtr_nirreps);
  moinfo.backtr_mo_active = init_int_array(moinfo.backtr_nirreps);
  moinfo.backtr_ao_first = init_int_array(moinfo.backtr_nirreps);
  moinfo.backtr_ao_last = init_int_array(moinfo.backtr_nirreps);
  moinfo.backtr_ao_orbspi = init_int_array(moinfo.backtr_nirreps);
  moinfo.backtr_ao_orbsym = init_int_array(moinfo.nao);

  moinfo.backtr_mo_first[0] = 0;
  moinfo.backtr_mo_last[0] = moinfo.nmo - 1;
  moinfo.backtr_mo_fstact[0] = 0;
  moinfo.backtr_mo_lstact[0] = moinfo.nmo - moinfo.nfzv - 1;
  moinfo.backtr_mo_orbspi[0] = moinfo.nmo;
  moinfo.backtr_mo_active[0] = moinfo.nmo - moinfo.nfzv;
  moinfo.backtr_ao_first[0] = 0;
  moinfo.backtr_ao_last[0] = moinfo.nao - 1;
  moinfo.backtr_ao_orbspi[0] = moinfo.nao;

  for (i=0; i<moinfo.nao; i++) moinfo.backtr_ao_orbsym[i] = 0;

}




void get_reorder_array(Options& options)
{
  int i, errcod;
  int **ras_opi;
  int j, k, l, fzv_offset;

  moinfo.order = init_int_array(moinfo.nmo);
  moinfo.order_alpha = init_int_array(moinfo.nmo);
  moinfo.order_beta = init_int_array(moinfo.nmo);

  /* for backtransforms, no reorder array...map Pitzer to Pitzer */
  if (params.wfn == "CI" || params.wfn == "DETCI" || params.wfn == "GVVPT2"
      || params.wfn == "MCSCF" || params.wfn == "OOCCD" || params.wfn == "ZAPTN"
      || params.wfn == "CASSCF" || params.wfn == "RASSCF" || params.wfn == "DETCAS") {

    ras_opi = init_int_matrix(MAX_RAS_SPACES,moinfo.nirreps);

    if (params.wfn == "GVVPT2" || params.wfn == "MCSCF")
      i=1;
    else i=0;

    if (!ras_set2(moinfo.nirreps, moinfo.nmo, params.fzc,
                 params.del_restr_docc, moinfo.orbspi,
                 moinfo.clsdpi, moinfo.openpi, moinfo.frdocc, moinfo.fruocc,
         moinfo.rstrdocc, moinfo.rstruocc, ras_opi, moinfo.order,
                 params.ras_type, i, options)) {
      fprintf(outfile, "Error in ras_set().  Aborting.\n");
      abort();
    }

    /* we just did this, but re-do it because ras_set() overwrites */
    if (params.treat_cor_as_fzc) {
      for (i=0; i<moinfo.nirreps; i++) {
        moinfo.frdocc[i] += moinfo.rstrdocc[i];
        moinfo.rstrdocc[i] = 0;
        if (params.backtr) {
          moinfo.fruocc[i] += moinfo.rstruocc[i];
          moinfo.rstruocc[i] = 0;
        }
      }
    }

    free_int_matrix(ras_opi);

  }

  else { /* default (CC, MP2, other) */
    reorder_qt(moinfo.clsdpi, moinfo.openpi, moinfo.frdocc, moinfo.fruocc,
           moinfo.order, moinfo.orbspi, moinfo.nirreps);
    reorder_qt_uhf(moinfo.clsdpi, moinfo.openpi, moinfo.frdocc, moinfo.fruocc,
           moinfo.order_alpha, moinfo.order_beta, moinfo.orbspi,
                   moinfo.nirreps);
  }

  /* Until I clean up all this, allow a PITZER flag */
  if(params.pitzer) {
    for(i=0; i < moinfo.nmo; i++){
      moinfo.order[i] = i;
      moinfo.order_alpha[i] = i;
      moinfo.order_beta[i] = i;
    }
  }


  if (params.print_reorder) {
    fprintf(outfile, "\nReordering array:\n");
    for (i=0; i<moinfo.nmo; i++) {
      fprintf(outfile, "%3d ", moinfo.order[i]);
    }
    fprintf(outfile, "\n");
  }

  /* construct an array to map the other direction, i.e., from correlated */
  /* to Pitzer order */
  moinfo.corr2pitz = init_int_array(moinfo.nmo);
  moinfo.corr2pitz_a = init_int_array(moinfo.nmo);
  moinfo.corr2pitz_b = init_int_array(moinfo.nmo);
  for (i=0; i<moinfo.nmo; i++) {
    j = moinfo.order[i];
    moinfo.corr2pitz[j] = i;

    j = moinfo.order_alpha[i];
    moinfo.corr2pitz_a[j] = i;

    j = moinfo.order_beta[i];
    moinfo.corr2pitz_b[j] = i;
  }

  if (params.backtr) {
    /* construct an array to map from correlated order to Pitzer order
     * less the frozen orbitals
     */
    fzv_offset = 0;

    moinfo.corr2pitz_nofzv = init_int_array(moinfo.nmo - moinfo.nfzv);
    moinfo.corr2pitz_nofzv_a = init_int_array(moinfo.nmo - moinfo.nfzv);
    moinfo.corr2pitz_nofzv_b = init_int_array(moinfo.nmo - moinfo.nfzv);
    for (i=0,k=0; i<moinfo.nirreps; i++) {
      for (j=0; j<moinfo.orbspi[i]; j++,k++) {
    if (j<moinfo.orbspi[i] - moinfo.fruocc[i]) {
      l = moinfo.order[k];
      moinfo.corr2pitz_nofzv[l] = k-fzv_offset;

      l = moinfo.order_alpha[k];
      moinfo.corr2pitz_nofzv_a[l] = k-fzv_offset;

      l = moinfo.order_beta[k];
      moinfo.corr2pitz_nofzv_b[l] = k-fzv_offset;
    }
    else fzv_offset++;
      }
    }

    if (params.print_reorder) {
      fprintf(outfile, "\nCorrelated to Pitzer (no fzv) Reordering Array:\n");
      for (i=0; i<moinfo.nmo-moinfo.nfzv; i++) {
    fprintf(outfile, "%3d ", moinfo.corr2pitz_nofzv[i]);
      }
      fprintf(outfile, "\n");
    }

  }
}



void get_one_electron_integrals()
{
  int i,stat;
  double enuc2,e_fzc;
  double *T, *V;

  if (params.backtr) return;

  moinfo.S = init_array(moinfo.noeints);
  T = init_array(moinfo.noeints);
  V = init_array(moinfo.noeints);
  moinfo.oe_ints = init_array(moinfo.noeints);
  if(params.ref == "UHF") {
    moinfo.fzc_operator_alpha = init_array(moinfo.noeints);
    moinfo.fzc_operator_beta = init_array(moinfo.noeints);
  }
  else moinfo.fzc_operator = init_array(moinfo.noeints);

  if (moinfo.S == NULL || T == NULL || V == NULL || moinfo.oe_ints == NULL) {
    printf("(transqt): Error mallocing one-electron ints\n");
    abort();
  }

  if (params.print_lvl)
    fprintf(outfile, "\n\tReading one-electron integrals...");
    stat = iwl_rdone(params.src_S_file,PSIF_SO_S,moinfo.S,moinfo.noeints,0,0,
                     outfile);
  if (!stat) {
    printf("(transqt): Error reading overlap ints\n");
    abort();
  }
  stat = iwl_rdone(params.src_T_file,PSIF_SO_T,T,moinfo.noeints,0,0,outfile);
  if (!stat) {
    printf("(transqt): Error reading kinetic energy ints\n");
    abort();
  }
  stat = iwl_rdone(params.src_V_file,PSIF_SO_V,V,moinfo.noeints,0,0,outfile);
  if (!stat) {
    printf("(transqt): Error reading potential energy ints\n");
    abort();
  }

  if (params.print_lvl) fprintf(outfile, "done.\n");

  for (i=0; i < moinfo.noeints; i++) {
    moinfo.oe_ints[i] = T[i] + V[i];

    if(params.ref == "UHF")
      moinfo.fzc_operator_alpha[i] = moinfo.fzc_operator_beta[i] = T[i] + V[i];
    else moinfo.fzc_operator[i] = T[i] + V[i];

  }

  free(T);
  free(V);

  /* if we're really doing frozen core, then we'll need the frozen core
   * density matrix when we're forming the two-electron contributions
   * to the frozen core operator, so go ahead and compute that now
   */
  if (params.fzc && moinfo.nfzc) {
    if(params.ref == "UHF") {
      moinfo.fzc_density_alpha = init_array(moinfo.noeints);
      fzc_density(moinfo.nirreps, moinfo.frdocc, moinfo.fzc_density_alpha,
          moinfo.scf_vector_alpha, moinfo.first, moinfo.first_so, moinfo.last_so, ioff);
      if (params.print_lvl > 2) {
    fprintf(outfile, "\nAlpha frozen core density matrix:\n");
    print_array(moinfo.fzc_density_alpha, moinfo.nso, outfile);
    fprintf(outfile, "\n");
      }
      moinfo.fzc_density_beta = init_array(moinfo.noeints);
      fzc_density(moinfo.nirreps, moinfo.frdocc, moinfo.fzc_density_beta,
          moinfo.scf_vector_beta, moinfo.first, moinfo.first_so, moinfo.last_so, ioff);
      if (params.print_lvl > 2) {
    fprintf(outfile, "\nBeta frozen core density matrix:\n");
    print_array(moinfo.fzc_density_beta, moinfo.nso, outfile);
    fprintf(outfile, "\n");
      }
    }
    else {
      moinfo.fzc_density = init_array(moinfo.noeints);
      if (params.ivo)
        ivo_density(moinfo.nirreps, moinfo.frdocc, moinfo.clsdpi,
                    moinfo.openpi, moinfo.fzc_density, moinfo.scf_vector,
                    moinfo.first, moinfo.first_so, moinfo.last_so, ioff);
      else
        fzc_density(moinfo.nirreps, moinfo.frdocc, moinfo.fzc_density,
                    moinfo.scf_vector, moinfo.first, moinfo.first_so,
                    moinfo.last_so, ioff);

      if (params.print_lvl > 2) {
    fprintf(outfile, "\nFrozen core density matrix:\n");
    print_array(moinfo.fzc_density, moinfo.nso, outfile);
    fprintf(outfile, "\n");
      }
    }
  }

}



/*
** construct_evects
**
** This function copies SCF eigenvectors from moinfo.scf_vector into
** an array of matrices.  Columns corresponding to inactive orbitals
** may be deleted if desired.
**
*/
double *** construct_evects(const char *spin, int nirreps, int *active, int *sopi,
                            int *orbspi, int *first_so, int *last_so,
                            int *first, int *last, int *fstact, int *lstact,
                            int printflag)
{

  int h, row, col, p, q;
  double ***evects;
  double **scf;

  std::string spin_str = spin;

  if(spin_str == "alpha")     scf = moinfo.scf_vector_alpha;
  else if(spin_str == "beta") scf = moinfo.scf_vector_beta;
  else if(spin_str == "RHF")  scf = moinfo.scf_vector;
  else throw PsiException ("Bad spin value in construct_evects", __FILE__, __LINE__);

  evects = (double ***) malloc(nirreps * sizeof(double **));

  for (h=0; h<nirreps; h++) {
    if (active[h]) {
      evects[h] = block_matrix(sopi[h],active[h]);
      row = -1;
      for(p=first_so[h]; p <= last_so[h]; p++) {
        row++; col = -1;
        for(q=fstact[h]; q <= lstact[h]; q++) {
          col++;
          evects[h][row][col] = scf[p][q];
        }
      }

      if(printflag) {
        fprintf(outfile,"\n\tMolecular Orbitals for Irrep %s\n",
                moinfo.labels[h]);
        print_mat(evects[h],sopi[h],active[h], outfile);
      }
    }

    else
      evects[h] = NULL;

  } /* end loop over irreps */

  return(evects);

}

void destruct_evects(int nirreps, double ***evects)
{
  int h;

  for (h=0; h<nirreps; h++)
    if (evects[h] != NULL) free_block(evects[h]);

  free(evects);
}


/*
** check_C
**
** This function checks the orthonormality of the C (SCF coefficient) matrix.
** We could do this one irrep at a time but it doesn't really matter.
*/
int check_C(int nso, int nmo, double **Cmat, double *S)
{
  double **Smat, **Tmat1, **Tmat2;
  double overlap_tol = 1.0E-12;
  int p, q, pq, failflag=0;

  Smat = block_matrix(nso, nso);
  Tmat1 = block_matrix(nso, nso);
  Tmat2 = block_matrix(nso, nso);

  for (p=0,pq=0; p<nso; p++)
    for (q=0; q<=p; q++,pq++)
      Smat[p][q] = Smat[q][p] = S[pq];

  /*
  fprintf(outfile, "C matrix:\n");
  print_mat(Cmat, nmo, nmo, outfile);
  fprintf(outfile, "S matrix:\n");
  print_mat(Smat, nmo, nmo, outfile);
  */

  /* T1 = C^+ S */
  mmult(Cmat, 1, Smat, 1, Tmat1, 0, nmo, nso, nso, 0);

  /* T2 = T1 C */
  mmult(Tmat1, 0, Cmat, 0, Tmat2, 0, nmo, nso, nmo, 0);

  /*
  fprintf(outfile, "C^+SC:\n");
  print_mat(Tmat2, nmo, nmo, outfile);
  */

  for (p=0; p<nmo; p++) {
    for (q=0; q<nmo; q++) {
      if (p==q) continue;
      if (fabs(Tmat2[p][q]) > overlap_tol) {
        fprintf(outfile, "(Check_C): C^+SC (%d,%d) = %lf\n", p, q, Tmat2[p][q]);
        failflag = 1;
      }
    }
  }

  for (p=0; p<nmo; p++) {
    if (fabs(Tmat2[p][p] - 1.00) > overlap_tol) {
      fprintf(outfile, "(Check_C): C^+SC (%d,%d) = %lf\n", p, p, Tmat2[p][p]);
      failflag = 1;
    }
  }

  free_block(Smat);
  free_block(Tmat1);
  free_block(Tmat2);

  return(failflag);
}

}} // end of namespace psi::transqt
