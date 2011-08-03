/*! \file
    \ingroup DETCI
    \brief Enter brief description of file here 
*/
/*
** PARAMS.CC: File contains functions which get or print the running
**    parameters for the CI calculation.
**
** David Sherrill, 16 November 1994
**
*/

#define EXTERN

#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <libipv1/ip_lib.h>
#include <libipv1/ip_data.gbl>
#include <libciomr/libciomr.h>
#include <libqt/qt.h>
#include <libchkpt/chkpt.h>
#include <psifiles.h>
#include "structs.h"
#include "globals.h"

namespace psi { namespace detci {

/*
** get_parameters(): Function gets the program running parameters such
**   as convergence.  These are stored in the Parameters data structure.
*/
void get_parameters(void)
{
   int i, j, k, errcod;
   int iopen=0, tval;
   char line1[133];
   double junk;
   
   /* default value of Ms0 depends on iopen but is modified below 
    * depending on value of opentype
    */

   /*
   chkpt_init(PSIO_OPEN_OLD);
   Parameters.Ms0 = !(chkpt_rd_iopen());
   chkpt_close();
   */

   /* need to figure out wheter to filter tei's */
   Parameters.dertype = options.get_str("DERTYPE");
   Parameters.wfn = options.get_str("WFN");

   /*
     Two-electron integrals: filter out what we don't need.  TRANSQT2
     supplies restricted orbitals always (well, for now).  It will also
     supply frozen core if it's a gradient calculation (need for orbital
     response) or an MCSCF (need for MO Hessian).  We normally want to
     filter all these out of the CI energy computation.  Likewise, we
     normally won't need restricted or frozen virtuals in the CI energy
     computation and should filter them out if they are in the TEI file
   */

  if (Parameters.dertype != "NONE" || Parameters.wfn == "DETCAS" ||
      Parameters.wfn == "CASSCF"   || Parameters.wfn == "RASSCF") {
    Parameters.filter_ints = 1;
  }
  else {
    Parameters.filter_ints = 0;
  }

 
  // CDS-TODO: Check these
  /* Parameters.print_lvl is set in detci.cc */
  /* Parameters.have_special_conv is set in detci.cc */

/* symbol * = put in read_options
 *        I = was internal parameter, not in read_options
 *        X = now it's an internal parameter, not in read_options
 */

  Parameters.ex_lvl = options.get_int("EX_LVL");
  Parameters.cc_ex_lvl = options.get_int("CC_EX_LVL");
  Parameters.val_ex_lvl = options.get_int("VAL_EX_LVL");
  Parameters.cc_val_ex_lvl = options.get_int("CC_VAL_EX_LVL");

  Parameters.cc_a_val_ex_lvl = -1;
  Parameters.cc_b_val_ex_lvl = -1;

  Parameters.num_roots = options.get_int("NUM_ROOTS");
 
  // CDS-TODO: Correct conversion of boolean to integer?
  Parameters.istop = options["ISTOP"].to_integer();
  Parameters.print_ciblks = options["PRINT_CIBLKS"].to_integer();

  Parameters.opentype = PARM_OPENTYPE_UNKNOWN;

  Parameters.ref_sym = options.get_int("REF_SYM");

  Parameters.oei_file = PSIF_OEI;  /* always need fzc operator */
  Parameters.tei_file = PSIF_MO_TEI;

  Parameters.hd_ave = EVANGELISTI;

  if (Parameters.wfn == "CASSCF")
    Parameters.fci = 1;
  else
    Parameters.fci = 0;

  Parameters.fci_strings = 0;

  if (Parameters.wfn == "ZAPTN") {
    Parameters.mpn = 1;
    Parameters.zaptn = 1;
  } 
  else {
    Parameters.mpn = 0;
    Parameters.zaptn = 0;
  }

  Parameters.z_scale_H = 0;
  Parameters.ras1_lvl = -1;
  Parameters.ras1_min = -1;
  Parameters.a_ras1_lvl = -1;
  Parameters.a_ras1_min = -1;
  Parameters.b_ras1_lvl = -1;
  Parameters.b_ras1_min = -1;
  Parameters.ras3_lvl = -1;
  Parameters.ras4_lvl = -1;

  if (Parameters.wfn == "DETCAS" ||
      Parameters.wfn == "CASSCF" ||
      Parameters.wfn == "RASSCF")
    Parameters.guess_vector = PARM_GUESS_VEC_DFILE;
  else
    Parameters.guess_vector = PARM_GUESS_VEC_H0_BLOCK;

  Parameters.update = UPDATE_DAVIDSON;

  Parameters.neg_only = 1;
  Parameters.nunits = 1;
  Parameters.first_tmp_unit = 50;
  Parameters.first_hd_tmp_unit = 0;
  Parameters.num_hd_tmp_units = 0;
  Parameters.first_c_tmp_unit = 0;
  Parameters.num_c_tmp_units = 0;
  Parameters.first_s_tmp_unit = 0;
  Parameters.num_s_tmp_units = 0;
  Parameters.first_d_tmp_unit = 0;
  Parameters.num_d_tmp_units = 0;
   
*  Parameters.restart = 0;
*  Parameters.bendazzoli = 0;

*  Parameters.cc = 0;
   errcod = ip_boolean("CC",&(Parameters.cc),0);

   if (Parameters.dertype == "FIRST" ||
      Parameters.wfn == "DETCAS" ||
      Parameters.wfn == "CASSCF" ||
      Parameters.wfn == "RASSCF")
   {
*    Parameters.convergence = 7;
*    Parameters.energy_convergence = 8;
*    Parameters.opdm = 1;
X    Parameters.opdm_write = 1;
*    Parameters.tpdm = 1;
X    Parameters.tpdm_write = 1;
*    Parameters.maxiter = 12;
   }

   else {
     if (Parameters.cc) {
       Parameters.convergence = 5;
       Parameters.energy_convergence = 7;
       Parameters.maxiter = 20;
     }
     else {
       Parameters.convergence = 4;
       Parameters.energy_convergence = 6;
       Parameters.maxiter = 12;
     }
     Parameters.opdm = 0;
     Parameters.opdm_write = 0;
     Parameters.tpdm = 0;
     Parameters.tpdm_write = 0;
   }

X  Parameters.opdm_file = PSIF_MO_OPDM;
*  Parameters.opdm_print = 0;
X  Parameters.opdm_diag = 0;
*  Parameters.opdm_wrtnos = 0;
X  Parameters.opdm_orbsfile = 76;
*  Parameters.opdm_ave = 0;
*  Parameters.opdm_orbs_root = -1;
*  Parameters.opdm_ke = 0;
X  Parameters.tpdm_file = PSIF_MO_TPDM;
*  Parameters.tpdm_print = 0;
*  Parameters.root = 1;

*  Parameters.nthreads = 1;
*  Parameters.sf_restrict = 0;
*  Parameters.print_sigma_overlap = 0;

   tval = 1;
   if (ip_exist("MULTP",0)) {
*    errcod = ip_data("MULTP","%d",&tval,0);
   }
*  Parameters.S = (((double) tval) - 1.0) / 2.0;
   errcod = ip_data("S","%lf",&(Parameters.S),0);

   errcod = ip_data("EX_LVL","%d",&(Parameters.ex_lvl),0);
   errcod = ip_data("CC_EX_LVL","%d",&(Parameters.cc_ex_lvl),0);
   errcod = ip_data("VAL_EX_LVL","%d",&(Parameters.val_ex_lvl),0);
   errcod = ip_data("CC_VAL_EX_LVL","%d",&(Parameters.cc_val_ex_lvl),0);
   errcod = ip_data("MAXITER","%d",&(Parameters.maxiter),0);
   errcod = ip_data("NUM_ROOTS","%d",&(Parameters.num_roots),0);
   errcod = ip_boolean("ISTOP",&(Parameters.istop),0);
   errcod = ip_data("PRINT","%d",&(Parameters.print_lvl),0);
   errcod = ip_boolean("PRINT_CIBLKS",&(Parameters.print_ciblks),0);
   errcod = ip_data("CONVERGENCE","%d",&(Parameters.convergence),0);
   errcod = ip_data("ENERGY_CONVERGENCE","%d",
               &(Parameters.energy_convergence),0);

   /* this handles backwards compatibility for PSI2 */
   errcod = ip_data("OPENTYPE","%s",line1,0);
   if (errcod == IPE_OK) {
      if (strcmp(line1, "NONE")==0) { 
         Parameters.opentype = PARM_OPENTYPE_NONE;
         Parameters.Ms0 = 1;
         }
      else if (strcmp(line1, "HIGHSPIN")==0) {
         Parameters.opentype = PARM_OPENTYPE_HIGHSPIN;
         Parameters.Ms0 = 0;
         }
      else if (strcmp(line1, "SINGLET")==0) {
         Parameters.opentype = PARM_OPENTYPE_SINGLET;
         Parameters.Ms0 = 1;
         }
      else Parameters.opentype = PARM_OPENTYPE_UNKNOWN;
      } 
   else { /* this handles new PSI3 keywords */
     errcod = ip_data("REFERENCE","%s",line1,0);
     /* fprintf(outfile, "line1 = "); 
     for (int ij = 0; ij<3; ij++) fprintf(outfile, "%1c", line1[ij]); 
     fprintf(outfile,"\n");
     */
     if (errcod == IPE_OK) {
       if (strcmp(line1, "RHF")==0) {
         Parameters.opentype = PARM_OPENTYPE_NONE;
         Parameters.Ms0 = 1;
         }
       else if (strcmp(line1, "ROHF")==0) {
         if (ip_data("MULTP","%d",&tval,0) == IPE_OK) {
           if (tval == 1) {
             Parameters.opentype = PARM_OPENTYPE_SINGLET;
             Parameters.Ms0 = 1;
             }
           else {
             Parameters.opentype = PARM_OPENTYPE_HIGHSPIN;
             Parameters.Ms0 = 0;
             }
           }
         else {
           fprintf(outfile, "detci: trouble reading MULTP\n");
           exit(0);
           }
         }  /* end ROHF parsing */
       else {
         fprintf(outfile, "detci: can only handle RHF or ROHF\n");
         exit(0);
         }
       }
     else Parameters.opentype = PARM_OPENTYPE_UNKNOWN;
     } /* end PSI3 parsing */

*  errcod = ip_boolean("MS0",&(Parameters.Ms0),0);
   errcod = ip_data("REF_SYM","%d",&(Parameters.ref_sym),0);
   errcod = ip_data("OEI_FILE","%d",&(Parameters.oei_file),0);
   errcod = ip_data("TEI_FILE","%d",&(Parameters.tei_file),0);
   errcod = ip_data("H0_BLOCKSIZE","%d",&(Parameters.h0blocksize),0);
   Parameters.h0guess_size = Parameters.h0blocksize;
   errcod = ip_data("H0_GUESS_SIZE","%d",&(Parameters.h0guess_size),0);
   if (Parameters.h0guess_size > Parameters.h0blocksize)
     Parameters.h0guess_size = Parameters.h0blocksize; 
   errcod = ip_data("H0_BLOCK_COUPLING_SIZE","%d",
                    &(Parameters.h0block_coupling_size),0);
   errcod = ip_boolean("H0_BLOCK_COUPLING",&(Parameters.h0block_coupling),0);
   errcod = ip_data("NPRINT","%d",&(Parameters.nprint),0);
   errcod = ip_data("CC_NPRINT","%d",&(Parameters.cc_nprint),0);
   errcod = ip_boolean("DETCI_FREEZE_CORE",&(Parameters.fzc),0);

   //CDS-TODO: Bool to Int
   if (options["FCI"].has_changed())
     Parameters.fci = options["FCI"].to_integer();

   if (Parameters.fci) Parameters.fci_strings = 1;
   errcod = ip_boolean("FCI_STRINGS",&(Parameters.fci_strings),0);
   errcod = ip_boolean("MIXED",&(Parameters.mixed),0);
   errcod = ip_boolean("MIXED4",&(Parameters.mixed4),0);
   errcod = ip_boolean("R4S",&(Parameters.r4s),0);
   errcod = ip_boolean("REPL_OTF",&(Parameters.repl_otf),0);
   errcod = ip_boolean("CALC_SSQ",&(Parameters.calc_ssq),0);
*  errcod = ip_boolean("MPN",&(Parameters.mpn),0);
   if (Parameters.mpn) {
     Parameters.fci = 1;
     Parameters.mpn_schmidt = FALSE;
     Parameters.wigner = TRUE;
     Parameters.guess_vector = PARM_GUESS_VEC_UNIT;
     Parameters.hd_ave = ORB_ENER;
     Parameters.update = UPDATE_DAVIDSON;
     Parameters.hd_otf = TRUE;
     Parameters.nodfile = TRUE;
     }
   errcod = ip_data("SAVE_MPN2","%d",&(Parameters.save_mpn2),0);
   errcod = ip_data("PERTURBATION_PARAMETER","%lf",
            &(Parameters.perturbation_parameter),0);
   
   if (Parameters.perturbation_parameter <= 1.0 && 
       Parameters.perturbation_parameter >= -1.0) Parameters.z_scale_H = 1;
/*
   else { fprintf(outfile, "Parameters.perturbation_parameters beyond the"
                 "bounds of -1.0 >= z <= 1.0\n");
         exit(0);
        }
*/
   errcod = ip_boolean("MPN_SCHMIDT",&(Parameters.mpn_schmidt),0);
   errcod = ip_boolean("WIGNER",&(Parameters.wigner),0);

   errcod = ip_data("A_RAS3_MAX","%d",&(Parameters.a_ras3_max),0);
   errcod = ip_data("B_RAS3_MAX","%d",&(Parameters.b_ras3_max),0);
   errcod = ip_data("RAS3_MAX","%d",&(Parameters.ras3_max),0);
   errcod = ip_data("RAS4_MAX","%d",&(Parameters.ras4_max),0);
   errcod = ip_data("RAS34_MAX","%d",&(Parameters.ras34_max),0);
   errcod = ip_data("CC_A_RAS3_MAX","%d",&(Parameters.cc_a_ras3_max),0);
   errcod = ip_data("CC_B_RAS3_MAX","%d",&(Parameters.cc_b_ras3_max),0);
   errcod = ip_data("CC_RAS3_MAX","%d",&(Parameters.cc_ras3_max),0);
   errcod = ip_data("CC_RAS4_MAX","%d",&(Parameters.cc_ras4_max),0);
   errcod = ip_data("CC_RAS34_MAX","%d",&(Parameters.cc_ras34_max),0);

   if (options["GUESS_VECTOR"].has_changed()) {
     std::string line1 = options.get_str("GUESS_VECTOR");
     if (line1 == "UNIT") 
       Parameters.guess_vector = PARM_GUESS_VEC_UNIT;
     else if (line1 == "H0_BLOCK") 
       Parameters.guess_vector = PARM_GUESS_VEC_H0_BLOCK;
     else if (line1 == "DFILE") 
       Parameters.guess_vector = PARM_GUESS_VEC_DFILE;
     /* else if (Parameters.mpn) Parameters.guess_vector = PARM_GUESS_VEC_UNIT; */ 
     else if (line1 == "IMPORT")
       Parameters.guess_vector = PARM_GUESS_VEC_IMPORT;
     else Parameters.guess_vector = PARM_GUESS_VEC_UNIT;
   }

   errcod = ip_data("ICORE", "%d", &(Parameters.icore),0);

   if (options["HD_AVE"].has_changed()) {
     std::string line1 = options.get_str("HD_AVE");
     if (line1 == "HD_EXACT")    Parameters.hd_ave = HD_EXACT;
     if (line1 == "HD_KAVE")     Parameters.hd_ave = HD_KAVE;
     if (line1 == "ORB_ENER")    Parameters.hd_ave = ORB_ENER;
     if (line1 == "EVANGELISTI") Parameters.hd_ave = EVANGELISTI;
     if (line1 == "LEININGER")   Parameters.hd_ave = LEININGER;
     if (line1 == "Z_KAVE")      Parameters.hd_ave = Z_HD_KAVE;
     /* if (Parameters.mpn) Parameters.hd_ave = ORB_ENER; */ 
   }

   errcod = ip_boolean("HD_OTF",&(Parameters.hd_otf),0); 
   if (errcod == IPE_OK) Parameters.hd_otf = TRUE;
   errcod = ip_boolean("NODFILE",&(Parameters.nodfile),0); 
   if (Parameters.num_roots > 1) Parameters.nodfile = FALSE;

   Parameters.diag_method = METHOD_DAVIDSON_LIU_SEM;
   if (options["DIAG_METHOD"].has_changed()) {
     std::string line1 = options.get_str("DIAG_METHOD");
     if (line1 == "RSP") Parameters.diag_method = METHOD_RSP;
     if (line1 == "OLSEN") Parameters.diag_method = METHOD_OLSEN;
     if (line1 == "MITRUSHENKOV") 
       Parameters.diag_method = METHOD_MITRUSHENKOV;
     if (line1 == "DAVIDSON") 
       Parameters.diag_method = METHOD_DAVIDSON_LIU_SEM;
     if (line1 == "SEM") 
       Parameters.diag_method = METHOD_DAVIDSON_LIU_SEM;
     if (line1 == "SEMTEST") 
       Parameters.diag_method = METHOD_RSPTEST_OF_SEM;
   }

   Parameters.precon = PRECON_DAVIDSON;
   if (options["PRECONDITIONER"].has_changed()) {
     std::string line1 = options.get_str("PRECONDITIONER");
     if (line1 == "LANCZOS") Parameters.precon = PRECON_LANCZOS;
     if (line1 == "DAVIDSON") Parameters.precon = PRECON_DAVIDSON;
     if (line1 == "GEN_DAVIDSON") 
       Parameters.precon = PRECON_GEN_DAVIDSON;
     if (line1 == "H0BLOCK") Parameters.precon = PRECON_GEN_DAVIDSON;
     if (line1 == "H0BLOCK_INV") 
       Parameters.precon = PRECON_H0BLOCK_INVERT;
     if (line1 == "ITER_INV") 
       Parameters.precon = PRECON_H0BLOCK_ITER_INVERT;
     if (line1 == "H0BLOCK_COUPLING") 
       Parameters.precon = PRECON_H0BLOCK_COUPLING;
     if (line1 == "EVANGELISTI") 
       Parameters.precon = PRECON_EVANGELISTI;
   }

   if (options["UPDATE"].has_changed()) {
     std::string line1 = options.get_str("UPDATE");
     if (line1 == "DAVIDSON") Parameters.update = UPDATE_DAVIDSON;
     if (line1 == "OLSEN") Parameters.update = UPDATE_OLSEN;
   }

   if (Parameters.diag_method < METHOD_DAVIDSON_LIU_SEM && 
       Parameters.update==UPDATE_DAVIDSON) {
     fprintf(outfile,"DAVIDSON update not available for OLSEN or MITRUSH"
             " iterators\n");
     Parameters.update = UPDATE_OLSEN;
     }
   if (Parameters.precon==PRECON_EVANGELISTI && (Parameters.update!=UPDATE_DAVIDSON 
        || Parameters.diag_method!=METHOD_DAVIDSON_LIU_SEM)) {
     fprintf(outfile,"EVANGELISTI preconditioner not available for OLSEN or"
                     " MITRUSH iterators or updates.\n");
     Parameters.update = UPDATE_DAVIDSON;
     }
   
   // errcod = ip_boolean("ZERO_BLOCKS",&(Parameters.zero_blocks),0);
   // if (Parameters.icore || !Parameters.mpn) Parameters.zero_blocks = 0;
   Parameters.num_init_vecs = Parameters.num_roots;
   errcod = ip_data("NUM_INIT_VECS","%d",&(Parameters.num_init_vecs),0);

   Parameters.collapse_size = options.get_int("COLLAPSE_SIZE");
   if (Parameters.collapse_size < 1) Parameters.collapse_size = 1;

   // CDS-TODO: Correct conversion of boolean to integer?
   Parameters.lse = options["LSE"].to_integer();

   Parameters.lse_collapse = options.get_int("LSE_COLLAPSE");
   if (Parameters.lse_collapse < 1) Parameters.lse_collapse = 3;

   Parameters.lse_tolerance = options.get_int("LSE_TOLERANCE");

   Parameters.maxnvect = options.get_int("MAXNVECT");

   if (Parameters.maxnvect == 0 &&  
       Parameters.diag_method == METHOD_DAVIDSON_LIU_SEM) {
      Parameters.maxnvect = Parameters.maxiter * Parameters.num_roots
         + Parameters.num_init_vecs;
      }
   else if (Parameters.maxnvect == 0 && 
            Parameters.diag_method == METHOD_RSPTEST_OF_SEM) {
      Parameters.maxnvect = Parameters.maxiter * Parameters.num_roots
         + Parameters.num_init_vecs;
      }
   else if (Parameters.maxnvect == 0 && 
            Parameters.diag_method == METHOD_MITRUSHENKOV) {
      Parameters.maxnvect = 2;
      }
   else if (Parameters.maxnvect == 0 && 
            Parameters.diag_method == METHOD_OLSEN) {
      Parameters.maxnvect = 1;
      }
   else { /* the user tried to specify a value for maxnvect...check it */
   /*    if (Parameters.maxnvect / (Parameters.collapse_size * 
         Parameters.num_roots) < 2) {
         fprintf(outfile, "maxnvect must be at least twice collapse_size *");
         fprintf(outfile, " num_roots.\n"); 
         exit(0);
         }
   */
      }
      
   errcod = ip_data("NUNITS", "%d", &(Parameters.nunits),0);
   errcod = ip_data("FIRST_TMP_UNIT", "%d", &(Parameters.first_tmp_unit),0);
   errcod = ip_data("FIRST_HD_TMP_UNIT","%d",&(Parameters.first_hd_tmp_unit),0);
   errcod = ip_data("FIRST_C_TMP_UNIT","%d",&(Parameters.first_c_tmp_unit),0);
   errcod = ip_data("FIRST_S_TMP_UNIT","%d",&(Parameters.first_s_tmp_unit),0);
   errcod = ip_data("FIRST_D_TMP_UNIT","%d",&(Parameters.first_d_tmp_unit),0);
   errcod = ip_data("NUM_HD_TMP_UNITS","%d",&(Parameters.num_hd_tmp_units),0);
   errcod = ip_data("NUM_C_TMP_UNITS","%d",&(Parameters.num_c_tmp_units),0);
   errcod = ip_data("NUM_S_TMP_UNITS","%d",&(Parameters.num_s_tmp_units),0);
   errcod = ip_data("NUM_D_TMP_UNITS","%d",&(Parameters.num_d_tmp_units),0);

   if (Parameters.first_hd_tmp_unit == 0) 
     Parameters.first_hd_tmp_unit = Parameters.first_tmp_unit;
/*   if ( (Parameters.num_hd_tmp_units == 0) && (!Parameters.hd_otf) ) */
   if (Parameters.num_hd_tmp_units == 0)
     Parameters.num_hd_tmp_units = 1;
   if (Parameters.first_c_tmp_unit == 0) Parameters.first_c_tmp_unit = 
      Parameters.first_hd_tmp_unit + Parameters.num_hd_tmp_units;
   if (Parameters.num_c_tmp_units == 0) Parameters.num_c_tmp_units = 
      Parameters.nunits;
   if (Parameters.first_s_tmp_unit == 0) Parameters.first_s_tmp_unit = 
      Parameters.first_c_tmp_unit + Parameters.num_c_tmp_units;
   if (Parameters.num_s_tmp_units == 0) Parameters.num_s_tmp_units = 
      Parameters.nunits;
   if (Parameters.first_d_tmp_unit == 0) Parameters.first_d_tmp_unit =
      Parameters.first_s_tmp_unit + Parameters.num_s_tmp_units;
 /*  if ( (Parameters.num_d_tmp_units == 0) && (!Parameters.nodfile) ) */
   if (Parameters.num_d_tmp_units == 0) 
     Parameters.num_d_tmp_units = 1;

   errcod = ip_boolean("RESTART",&(Parameters.restart),0);
 /* obsolete due to new restart procedure
   errcod = ip_data("RESTART_ITER","%d",&(Parameters.restart_iter),0);
   errcod = ip_data("RESTART_VECS","%d",&(Parameters.restart_vecs),0);
   if (Parameters.restart && (errcod!=IPE_OK || Parameters.restart_vecs==0)) {
      fprintf(outfile, "For RESTART must specify nonzero RESTART_VECS\n");
      exit(0);
      }
 */ 
   errcod = ip_boolean("BENDAZZOLI",&(Parameters.bendazzoli),0) ;
   if (Parameters.bendazzoli & !Parameters.fci) Parameters.bendazzoli=0;

   /* Parse the OPDM stuff.  It is possible to give incompatible options,
    * but we will try to eliminate some of those.  Parameters_opdm will
    * function as the master switch for all other OPDM parameters.
    */
   errcod = ip_boolean("OPDM_PRINT",&(Parameters.opdm_print),0);
   // Make this an internal parameter
   // errcod = ip_data("OPDM_FILE","%d",&(Parameters.opdm_file),0);
   errcod = ip_boolean("WRTNOS",&(Parameters.opdm_wrtnos),0);
   // Make this an internal parameter, essentially same as WRTNOS
   // errcod = ip_boolean("OPDM_DIAG",&(Parameters.opdm_diag),0);
   errcod = ip_boolean("OPDM_AVE",&(Parameters.opdm_ave),0);
   // Make an internal parameter
   // errcod = ip_data("ORBSFILE","%d",&(Parameters.opdm_orbsfile),0);
   errcod = ip_data("ORBS_ROOT","%d",&(Parameters.opdm_orbs_root),0);
   errcod = ip_boolean("OPDM_KE",&(Parameters.opdm_ke),0);
   
   if (Parameters.opdm_orbs_root != -1) Parameters.opdm_orbs_root -= 1;
   if (Parameters.opdm_orbs_root < 0) Parameters.opdm_orbs_root = 0;
   if (Parameters.opdm_wrtnos) Parameters.opdm_diag = 1;
   if (Parameters.opdm_print || Parameters.opdm_diag || Parameters.opdm_wrtnos 
       || Parameters.opdm_ave || Parameters.opdm_ke) Parameters.opdm = 1;
   errcod = ip_boolean("OPDM",&(Parameters.opdm),0);
   if (Parameters.opdm) Parameters.opdm_write = 1;
//   No reason why the user would need to change this, make internal param
//   errcod = ip_boolean("OPDM_WRITE",&(Parameters.opdm_write),0);
   errcod = ip_boolean("OPDM_PRINT",&(Parameters.opdm_print),0);
   errcod = ip_data("OPDM_DIAG","%d",&(Parameters.opdm_diag),0);

   /* transition density matrices */
   Parameters.tdm_print = 0; Parameters.tdm_write = 0;
   Parameters.transdens = 0;
   errcod = ip_boolean("TDM_PRINT",&(Parameters.tdm_print),0);
   errcod = ip_boolean("TDM_WRITE",&(Parameters.tdm_write),0);
   if (Parameters.tdm_print || Parameters.tdm_write || 
       (Parameters.num_roots > 1))
*    Parameters.transdens = 1;
   else 
     Parameters.transdens = 0;
   errcod = ip_boolean("TRANSITION_DENSITY",&(Parameters.transdens),0);
   if (Parameters.transdens) Parameters.tdm_write = 1;
   errcod = ip_boolean("TDM_WRITE",&(Parameters.tdm_write),0);
  
   /* dipole or transition dipole moment? */
   if (Parameters.opdm) Parameters.dipmom = 1;
   else Parameters.dipmom = 0;
   if (Parameters.wfn == "RASSCF" || 
       Parameters.wfn == "CASSCF" ||
       Parameters.wfn == "DETCAS")
   Parameters.dipmom = 0;
   if (Parameters.transdens) Parameters.dipmom = 1;
*  errcod = ip_boolean("DIPMOM",&(Parameters.dipmom),0);
   if (Parameters.dipmom == 1) Parameters.opdm = 1; 
 
   errcod = ip_data("ROOT","%d",&(Parameters.root),0);
   Parameters.root -= 1;

   errcod = ip_boolean("TPDM",&(Parameters.tpdm),0);
   if (Parameters.tpdm) Parameters.tpdm_write = 1;
//   Made this an internal parameter
//   errcod = ip_boolean("TPDM_WRITE",&(Parameters.tpdm_write),0);
   errcod = ip_boolean("TPDM_PRINT",&(Parameters.tpdm_print),0);
//   Made this an internal parameter
//   errcod = ip_data("TPDM_FILE","%d",&(Parameters.tpdm_file),0);

   if (Parameters.guess_vector == PARM_GUESS_VEC_DFILE &&
       Parameters.wfn != "DETCAS" &&
       Parameters.wfn != "CASSCF" &&
       Parameters.wfn != "RASSCF")
     {

      chkpt_init(PSIO_OPEN_OLD);
      i = chkpt_rd_phase_check();
      chkpt_close();

      if (!i) {
         fprintf(outfile, "Can't use d file guess: SCF phase not checked\n");
         if (Parameters.h0guess_size) {
            Parameters.guess_vector = PARM_GUESS_VEC_H0_BLOCK;
            if (Parameters.precon == PRECON_GEN_DAVIDSON)
              Parameters.precon = PRECON_H0BLOCK_ITER_INVERT;
         }
         else Parameters.guess_vector = PARM_GUESS_VEC_UNIT;
      }
     }
   if (Parameters.num_init_vecs < Parameters.num_roots)
     Parameters.num_init_vecs = Parameters.num_roots;
   if (Parameters.guess_vector == PARM_GUESS_VEC_UNIT &&
       Parameters.num_init_vecs > 1) {
     Parameters.guess_vector = PARM_GUESS_VEC_H0_BLOCK;
     fprintf(outfile,"Warning: Unit vec option not available for more than"
             " one root\n");
     }
   if (Parameters.guess_vector == PARM_GUESS_VEC_UNIT)
     Parameters.h0blocksize = Parameters.h0guess_size = 1;

   errcod = ip_data("NTHREADS", "%d", &(Parameters.nthreads),0);
   if (Parameters.nthreads < 1) Parameters.nthreads = 1;

   Parameters.export_ci_vector = 0;
   errcod = ip_boolean("EXPORT_VECTOR", &(Parameters.export_ci_vector), 0);

   Parameters.num_export = 0;
   if (Parameters.export_ci_vector) {
     Parameters.num_export = 1;
     errcod = ip_data("NUM_EXPORT", "%d", &(Parameters.num_export), 0);
     if (Parameters.num_export > Parameters.num_roots) {
       fprintf(outfile, "Warning: can't export %d roots if %d requested\n",
         Parameters.num_export, Parameters.num_roots);
       Parameters.num_export = Parameters.num_roots;
     }
   } 
   
   errcod = ip_boolean("SF_RESTRICT",&(Parameters.sf_restrict),0);
   errcod = ip_boolean("SIGMA_OVERLAP",&(Parameters.print_sigma_overlap),0);
   
   if (Parameters.cc) Parameters.ex_lvl = Parameters.cc_ex_lvl + 2;

   Parameters.ex_allow = (int *)malloc(Parameters.ex_lvl*sizeof(int));
   if (ip_exist("EX_ALLOW",0)) {
     ip_count("EX_ALLOW", &i, 0);
     if (i != Parameters.ex_lvl) {
       fprintf(outfile,"Dim. of EX_ALLOW must be %d\n", 
               Parameters.ex_lvl);
       exit(0);
     }
     for (i=0;i<Parameters.ex_lvl;i++) {
*      errcod = ip_data("EX_ALLOW","%d",&(Parameters.ex_allow[i]),1,i);
     }
   }
   else {
     for (i=0;i<Parameters.ex_lvl;i++) {
       Parameters.ex_allow[i] = 1;
     }
   }

   /* The filter_guess options are used to filter out some trial
      vectors which may not have the appropriate phase convention
      between two determinants.  This is useful to remove, e.g.,
      delta states when a sigma state is desired.  The user
      inputs two determinants (by giving the absolute alpha string
      number and beta string number for each), and also the
      desired phase between these two determinants for guesses
      which are to be kept.
    */
   Parameters.filter_guess = 0;
*  errcod = ip_boolean("FILTER_GUESS",&(Parameters.filter_guess),0);
   if (errcod == IPE_OK && Parameters.filter_guess == 1) {
     Parameters.filter_guess_sign = 1;
*    errcod = ip_data("FILTER_GUESS_SIGN","%d",
                      &(Parameters.filter_guess_sign),0);
     if (errcod != IPE_OK || (Parameters.filter_guess_sign != 1 &&
         Parameters.filter_guess_sign != -1)) {
       fprintf(outfile, "FILTER_GUESS_SIGN should be 1 or -1 !\n");
       abort();
     }
*    errcod = ip_count("FILTER_GUESS_DET1",&i,0); 
     if (errcod != IPE_OK || i != 2) {
       fprintf(outfile, "Need to specify FILTER_GUESS_DET1 = "
                        "(alphastr betastr)\n");
       abort();
     }
     else {
       errcod = ip_data("FILTER_GUESS_DET1","%d",
                        &(Parameters.filter_guess_Ia),1,0);
       errcod = ip_data("FILTER_GUESS_DET1","%d",
                        &(Parameters.filter_guess_Ib),1,1);
     }

*    errcod = ip_count("FILTER_GUESS_DET2",&i,0); 
     if (errcod != IPE_OK || i != 2) {
       fprintf(outfile, "Need to specify FILTER_GUESS_DET2 = "
                        "(alphastr betastr)\n");
       abort();
     }
     else {
       errcod = ip_data("FILTER_GUESS_DET2","%d",
                        &(Parameters.filter_guess_Ja),1,0);
       errcod = ip_data("FILTER_GUESS_DET2","%d",
                        &(Parameters.filter_guess_Jb),1,1);
     }
   } /* end the filter_guess stuff */

   /* sometimes the filter guess stuff is not sufficient in that
      some states come in that we don't want.  We can help exclude
      them by explicitly zeroing out certain determinants, if that
      is correct for the desired state.  This stuff will allow the
      user to select a determinant which should always have a zero
      coefficient in the desired target state
    */
   Parameters.filter_zero_det = 0;
*  if (ip_exist("FILTER_ZERO_DET",0)) {
     errcod = ip_count("FILTER_ZERO_DET",&i,0); 
     if (errcod != IPE_OK || i != 2) {
       fprintf(outfile, "Need to specify FILTER_ZERO_DET = "
                        "(alphastr betastr)\n");
       abort();
     }
     Parameters.filter_zero_det = 1;
     errcod = ip_data("FILTER_ZERO_DET","%d",
                      &(Parameters.filter_zero_det_Ia),1,0);
     errcod = ip_data("FILTER_ZERO_DET","%d",
                      &(Parameters.filter_zero_det_Ib),1,1);
   }

   /* Does the user request a state-averaged calculation? */
*  if (ip_exist("AVERAGE_STATES",0)) {
     ip_count("AVERAGE_STATES", &i, 0);
     if (i < 1 || i > Parameters.num_roots) {
       fprintf(outfile,"Invalid number of states to average (%d)\n", i);
       exit(1);
     }
     Parameters.average_states = init_int_array(i);
     Parameters.average_weights = init_array(i);
     Parameters.average_num = i;
     for (i=0;i<Parameters.average_num;i++) {
       errcod = ip_data("AVERAGE_STATES","%d",
         &(Parameters.average_states[i]),1,i);
       if (Parameters.average_states[i] < 1) {
         fprintf(outfile,"AVERAGE_STATES start numbering from 1.\n");
         fprintf(outfile,"Invalid state number %d\n", 
           Parameters.average_states[i]);
         exit(1);
       }
       Parameters.average_states[i] -= 1; /* number from 1 externally */
       Parameters.average_weights[i] = 1.0/((double)Parameters.average_num);
     }
*    if (ip_exist("AVERAGE_WEIGHTS",0)) {
       ip_count("AVERAGE_WEIGHTS", &i, 0);
       if (i != Parameters.average_num) {
         fprintf(outfile,"Mismatched number of average weights (%d)\n", i);
         exit(0);
       }
       for (i=0; i<Parameters.average_num; i++) {
         errcod = ip_data("AVERAGE_WEIGHTS","%lf", 
           &(Parameters.average_weights[i]),1,i);
       }
     }
      
   if (Parameters.average_num > 1) Parameters.opdm_ave = 1;

     if ((!ip_exist("ROOT",0)) && (Parameters.average_num==1)) {
       Parameters.root = Parameters.average_states[0];
     }
   }

   else {
     Parameters.average_num = 1;
     Parameters.average_states = init_int_array(1);
     Parameters.average_weights = init_array(1);
     Parameters.average_states[0] = Parameters.root;
     Parameters.average_weights[0] = 1.0;
   } /* end state-average parsing */

   /* Follow a vector to determine the root number? */
   Parameters.follow_vec_num = 0;

   /* CDS-TODO: I'm going to rearrange this into 
      FOLLOW_VECTOR_ALPHAS, FOLLOW_VECTOR_BETAS, and FOLLOW_VECTOR_COEFS
   */
   if (ip_exist("FOLLOW_VECTOR",0)) {
     errcod = ip_count("FOLLOW_VECTOR",&i,0); 
     if (errcod != IPE_OK || (i % 2)!=0) {
       fprintf(outfile, "Need to specify FOLLOW_VECTOR = "
                        "((alphastr_i betastr_i) coeff_i ... )\n");
       abort();
     }
     i = i/2;
     Parameters.follow_vec_num = i;     
*    Parameters.follow_vec_coef   = init_array(i);
*    Parameters.follow_vec_Ia     = init_int_array(i);
*    Parameters.follow_vec_Ib     = init_int_array(i);
     Parameters.follow_vec_Iac    = init_int_array(i);
     Parameters.follow_vec_Ibc    = init_int_array(i);
     Parameters.follow_vec_Iaridx = init_int_array(i);
     Parameters.follow_vec_Ibridx = init_int_array(i);

     /* now parse each piece */
     for (i=0; i<Parameters.follow_vec_num; i++) {
       errcod = ip_data("FOLLOW_VECTOR","%d",
                        Parameters.follow_vec_Ia+i,2,i*2,0);
       if (errcod != IPE_OK) {
         fprintf(outfile, "Trouble parsing FOLLOW_VECTOR\n");
         abort();
       }
       errcod = ip_data("FOLLOW_VECTOR","%d",
                        Parameters.follow_vec_Ib+i,2,i*2,1);
       if (errcod != IPE_OK) {
         fprintf(outfile, "Trouble parsing FOLLOW_VECTOR\n");
         abort();
       }
       errcod = ip_data("FOLLOW_VECTOR","%lf",
                        Parameters.follow_vec_coef+i,1,i*2+1);
       if (errcod != IPE_OK) {
         fprintf(outfile, "Trouble parsing FOLLOW_VECTOR\n");
         abort();
       }
     } /* end loop over parsing */
   } /* end follow vector stuff */

   /* make sure SA weights add up to 1.0 */
   for (i=0,junk=0.0; i<Parameters.average_num; i++) {
     junk += Parameters.average_weights[i]; 
   }
   if (junk <= 0.0) {
     fprintf(outfile, "Error: AVERAGE WEIGHTS add up to %12.6lf\n", junk);
     exit(0);
   }
   for (i=0; i<Parameters.average_num; i++) {
     Parameters.average_weights[i] /= junk; 
   }

  Parameters.cc_export = 0;
  Parameters.cc_import = 0;
* errcod = ip_boolean("CC_EXPORT", &(Parameters.cc_export), 0);
* errcod = ip_boolean("CC_IMPORT", &(Parameters.cc_import), 0);

  Parameters.cc_fix_external = 0;
* errcod = ip_boolean("CC_FIX_EXTERNAL", &(Parameters.cc_fix_external), 0);
  Parameters.cc_fix_external_min = 1;
* errcod = ip_data("CC_FIX_EXTERNAL_MIN", "%d",
    &(Parameters.cc_fix_external_min),0);

  Parameters.cc_variational = 0;
* errcod = ip_boolean("CC_VARIATIONAL", &(Parameters.cc_variational), 0);

  Parameters.cc_mixed = 1;
* errcod = ip_boolean("CC_MIXED", &(Parameters.cc_mixed), 0);

  /* update using orb eigvals or not? */
  Parameters.cc_update_eps = 1;
* errcod = ip_boolean("CC_UPDATE_EPS", &(Parameters.cc_update_eps), 0);

* Parameters.diis = 1;             /* only kicks in for CC anyway */
* Parameters.diis_start = 1;       /* iteration to turn on DIIS */
* Parameters.diis_freq  = 1;       /* how often to do a DIIS extrapolation */
* Parameters.diis_min_vecs = 2;
* Parameters.diis_max_vecs = 5;

  errcod = ip_boolean("DIIS", &(Parameters.diis), 0);
  errcod = ip_data("DIIS_START","%d",&(Parameters.diis_start),0);
  errcod = ip_data("DIIS_FREQ","%d",&(Parameters.diis_freq),0);
  errcod = ip_data("DIIS_MIN_VECS","%d",&(Parameters.diis_min_vecs),0);
  errcod = ip_data("DIIS_MAX_VECS","%d",&(Parameters.diis_max_vecs),0);
   
  /* parse cc_macro = (
       (ex_lvl max_holes_I max_parts_IV max_I+IV)
       ...
     )
     This says to eliminate T blocks in which: [the excitation level
     (holes in I + II) is equal to ex_lvl] AND [there are more than
     max_holes_I holes in RAS I, there are more than max_parts_IV
     particles in RAS IV, OR there are more than max_I+IV quasiparticles
     in RAS I + RAS IV]
   */
  // CDS-TODO: Here we need an array of arrays.  This seems like it
  // should work by analogy to src/bin/cis/local.cc.  Check.
  // Here's the CIS example:
  /*
  if (options["DOMAINS"].size() > 0) {
    num_entries = options["DOMAINS"].size();
    for(i=0; i<num_entries; i++) {
      entry_len = options["DOMAINS"][i].size();
      orbital = options["DOMAINS"][i][0].to_integer();

      for(j=1; j<entry_len; j++) {
        atom = options["DOMAINS"][i][j].to_integer();
        domain[orbital][atom] = 1;
        domain_len[orbital]++;
      }
    }
  }
  */

  Parameters.cc_macro_on = 0;
* if (Parameters.cc && ip_exist("CC_MACRO",0)) {
    errcod = ip_count("CC_MACRO",&i,0);
    Parameters.cc_macro = init_int_matrix(Parameters.cc_ex_lvl +
      Parameters.cc_val_ex_lvl + 1, 3);
    Parameters.cc_macro_parsed = init_int_array(Parameters.cc_ex_lvl +
      Parameters.cc_val_ex_lvl + 1);
    for (j=0; j<i; j++) {
      errcod = ip_data("CC_MACRO","%d",&k,2,j,0);
      if (errcod != IPE_OK) {
        fprintf(outfile, "Trouble parsing CC_MACRO\n");
        abort();
      }
      errcod = ip_data("CC_MACRO","%d",Parameters.cc_macro[k],2,j,1);
      if (errcod != IPE_OK) {
        fprintf(outfile, "Trouble parsing CC_MACRO\n");
        abort();
      }
      errcod = ip_data("CC_MACRO","%d",Parameters.cc_macro[k]+1,2,j,2);
      if (errcod != IPE_OK) {
        fprintf(outfile, "Trouble parsing CC_MACRO\n");
        abort();
      }
      errcod = ip_data("CC_MACRO","%d",Parameters.cc_macro[k]+2,2,j,3);
      if (errcod != IPE_OK) {
        fprintf(outfile, "Trouble parsing CC_MACRO\n");
        abort();
      }
      Parameters.cc_macro_parsed[k] = 1;
    }
    Parameters.cc_macro_on = 1;
  } /* end parsing of CC_MACRO */

}


/*
** print_parameters(): Function prints the program's running parameters
**   found in the Parameters structure.
*/
void print_parameters(void)
{
   int i;

   fprintf(outfile, "\n");
   fprintf(outfile, "PARAMETERS: \n");
   fprintf(outfile, "   EX LVL        =   %6d      H0 BLOCKSIZE =   %6d\n", 
      Parameters.ex_lvl, Parameters.h0blocksize);
   fprintf(outfile, "   VAL EX LVL    =   %6d      H0 GUESS SIZE=   %6d\n", 
      Parameters.val_ex_lvl, Parameters.h0guess_size);
   fprintf(outfile, "   H0COUPLINGSIZE=   %6d      H0 COUPLING  =   %6s\n", 
      Parameters.h0block_coupling_size, Parameters.h0block_coupling ? "yes" : "no");
   fprintf(outfile, "   NPRINT        =   %6d\n", Parameters.nprint);
   fprintf(outfile, "   MAXITER       =   %6d      FREEZE CORE  =   %6s\n", 
      Parameters.maxiter, Parameters.fzc ? "yes" : "no");
   fprintf(outfile, "   NUM ROOTS     =   %6d      ICORE        =   %6d\n", 
      Parameters.num_roots, Parameters.icore);
   fprintf(outfile, "   PRINT         =   %6d      FCI          =   %6s\n", 
      Parameters.print_lvl, Parameters.fci ? "yes" : "no");
   if (Parameters.have_special_conv) 
      fprintf(outfile, 
         "   CONV          =   %8.2g    MIXED        =   %6s\n", 
         Parameters.special_conv, Parameters.mixed ? "yes" : "no");
   else
      fprintf(outfile, "   CONV          =   %6d      MIXED        =   %6s\n", 
         Parameters.convergence, Parameters.mixed ? "yes" : "no");

   fprintf(outfile, "   E CONV        =   %6d      MIXED4       =   %6s\n", 
      Parameters.energy_convergence, Parameters.mixed4 ? "yes" : "no");
   fprintf(outfile, "   OEI FILE      =   %6d      R4S          =   %6s\n", 
      Parameters.oei_file, Parameters.r4s ? "yes" : "no");
   fprintf(outfile, "   REPL OTF     =   %6s\n",  
      Parameters.repl_otf ? "yes" : "no");
   fprintf(outfile, "   TEI FILE      =   %6d      DIAG METHOD  =   ", 
      Parameters.tei_file);

   switch (Parameters.diag_method) {
      case 0:
         fprintf(outfile, "%6s\n", "RSP");
         break;
      case 1:
         fprintf(outfile, "%6s\n", "OLSEN");
         break;
      case 2:
         fprintf(outfile, "%6s\n", "MITRUS");
         break;
      case 3:
         fprintf(outfile, "%6s\n", "SEM");
         break;
      case 4:
         fprintf(outfile, "%6s\n", "SEMTEST");
         break;
      default:
         fprintf(outfile, "%6s\n", "???");
         break;
      } 

   fprintf(outfile, "   PRECONDITIONER= ");
   switch (Parameters.precon) {
      case PRECON_LANCZOS:
         fprintf(outfile, "%6s", " LANCZOS    ");
         break;
      case PRECON_DAVIDSON:
         fprintf(outfile, "%6s", "DAVIDSON    ");
         break;
      case PRECON_GEN_DAVIDSON:
         fprintf(outfile, "%6s", "GEN_DAVIDSON");
         break;
      case PRECON_H0BLOCK_INVERT:
         fprintf(outfile, "%6s", "H0BLOCK_INV ");
         break;
      case PRECON_H0BLOCK_ITER_INVERT:
         fprintf(outfile, "%6s", "ITER_INV    ");
         break;
      case PRECON_H0BLOCK_COUPLING:
         fprintf(outfile, "%6s", "H0_COUPLING ");
         break;
      case PRECON_EVANGELISTI:
         fprintf(outfile, "%6s", "EVANGELISTI ");
         break;
      default:
         fprintf(outfile, "%6s", "???         ");
         break;
      } 

   fprintf(outfile, "  UPDATE       =   ");
   switch (Parameters.update) {
     case 1:
       fprintf(outfile, "%6s\n", "DAVIDSON");
       break;
     case 2:
       fprintf(outfile, "%6s\n", "OLSEN");  
       break;
     default:
       fprintf(outfile, "%6s\n", "???");
       break;
      }

   fprintf(outfile, "   S             =   %6lf     Ms0          =   %6s\n",
      Parameters.S, Parameters.Ms0 ? "yes" : "no");           
   fprintf(outfile, "   MAXNVECT     =   %6d\n", Parameters.maxnvect);
   fprintf(outfile, "   RESTART       =   %6s\n",
      Parameters.restart ? "yes" : "no");
   fprintf(outfile, "   GUESS VECTOR  =  ");
   switch (Parameters.guess_vector) {
      case PARM_GUESS_VEC_UNIT:
         fprintf(outfile, "%7s", "UNIT");
         break; 
      case PARM_GUESS_VEC_H0_BLOCK:
         fprintf(outfile, "%7s", "H0BLOCK");
         break;
      case PARM_GUESS_VEC_DFILE:
         fprintf(outfile, "%7s", "D FILE");
         break;
      case PARM_GUESS_VEC_IMPORT:
         fprintf(outfile, "%7s", "IMPORT");
         break;
      default:
         fprintf(outfile, "%7s", "???");
         break;
      }
   fprintf(outfile, "      OPENTYPE     = ");
   switch (Parameters.opentype) {
      case PARM_OPENTYPE_NONE:
         fprintf(outfile, "%8s\n", "NONE");
         break;
      case PARM_OPENTYPE_HIGHSPIN:
         fprintf(outfile, "%8s\n", "HIGHSPIN");
         break;
      case PARM_OPENTYPE_SINGLET:
         fprintf(outfile, "%8s\n", "SINGLET");
         break;
      default:
         fprintf(outfile, "%8s\n", "???"); 
         break;
      }
   if (Parameters.ref_sym == -1)
      fprintf(outfile, "      REF SYM      =   %6s\n", "auto");
   else
      fprintf(outfile, "      REF SYM      =   %6d\n", Parameters.ref_sym);

   fprintf(outfile, "   COLLAPSE SIZE =   %6d", Parameters.collapse_size);
   fprintf(outfile, "      HD AVE       =");
   switch (Parameters.hd_ave) {
     case HD_EXACT:
       fprintf(outfile," %11s\n", "HD_EXACT");
       break;
     case HD_KAVE:
       fprintf(outfile," %11s\n", "HD_KAVE");
       break;
     case ORB_ENER:
       fprintf(outfile," %11s\n", "ORB_ENER");
       break;
     case EVANGELISTI:
       fprintf(outfile," %11s\n", "EVANGELISTI");
       break;
     case LEININGER:
       fprintf(outfile," %11s\n", "LEININGER");
       break;
     default:
       fprintf(outfile," %11s\n", "???");       
       break;
     }

   fprintf(outfile, "   LSE           =   %6s      LSE ITER     =   %6d\n", 
           Parameters.lse ? "yes" : "no", Parameters.lse_iter);
   fprintf(outfile, "   HD OTF        =   %6s      NO DFILE     =   %6s\n", 
           Parameters.hd_otf ? "yes" : "no", Parameters.nodfile ? "yes":"no");
   fprintf(outfile, "   MPN           =   %6s      MPN SCHMIDT  =   %6s\n",
           Parameters.mpn ? "yes":"no", Parameters.mpn_schmidt ? "yes":"no");
   fprintf(outfile, " ZAPTN           =   %6s      WIGNER       =   %6s\n",
           Parameters.zaptn ? "yes":"no", Parameters.wigner ? "yes":"no");
   fprintf(outfile, "   PERT Z        =   %1.4f      ROOT         =   %6d\n",
           Parameters.perturbation_parameter, Parameters.root);
   fprintf(outfile, "   NTHREADS     =   %6d\n",
           Parameters.nthreads);
   fprintf(outfile, "   EXPORT VECTOR =   %6s      NUM EXPORT   =   %6d\n",
           Parameters.export_ci_vector ? "yes":"no", Parameters.num_export);
   fprintf(outfile, "   FILTER_GUESS  =   %6s      SF_RESTRICT  =   %6s\n",
           Parameters.filter_guess ?  "yes":"no",
	   Parameters.sf_restrict ? "yes":"no");
   if (Parameters.cc && Parameters.diis) {
     fprintf(outfile, "   DIIS START    =   %6d      DIIS FREQ    =   %6d\n",
             Parameters.diis_start, Parameters.diis_freq);
     fprintf(outfile, "   DIIS MIN VECS =   %6d      DIIS MAX VECS=   %6d\n",
             Parameters.diis_min_vecs, Parameters.diis_max_vecs);
   }
   fprintf(outfile, "   OPDM          =   %6s      TRANS DENSITY=   %6s\n",
           Parameters.opdm ?  "yes":"no",
	   Parameters.transdens ? "yes":"no");
   fprintf(outfile, "\n   FILES         = %3d %2d %2d %2d\n",
      Parameters.first_hd_tmp_unit, Parameters.first_c_tmp_unit,
      Parameters.first_s_tmp_unit, Parameters.first_d_tmp_unit);
   
   fprintf(outfile, "\n   EX_ALLOW      = ");
   for (i=0;i<Parameters.ex_lvl;i++) {
     fprintf(outfile, "%2d ", Parameters.ex_allow[i]);
   }

   fprintf(outfile, "\n   STATE AVERAGE = ");
   for (i=0; i<Parameters.average_num; i++) {
     if (i%5==0 && i!=0) fprintf(outfile, "\n");
     fprintf(outfile, "%2d(%4.2lf) ",Parameters.average_states[i]+1,
       Parameters.average_weights[i]);
   }

   fprintf(outfile, "\n   STATE AVERAGE = ");
   for (i=0; i<Parameters.average_num; i++) {
     if (i%5==0 && i!=0) fprintf(outfile, "\n");
     fprintf(outfile, "%2d(%4.2lf) ",Parameters.average_states[i]+1,
       Parameters.average_weights[i]);
   }

   if (Parameters.follow_vec_num > 0) {
     fprintf(outfile,"\nDensity matrices will follow vector like:\n");
     for (i=0; i<Parameters.follow_vec_num; i++) 
       fprintf(outfile, "(%d %d) %12.6lf\n", Parameters.follow_vec_Ia[i],
         Parameters.follow_vec_Ib[i], Parameters.follow_vec_coef[i]);
   }

   fprintf(outfile, "\n\n");
   fflush(outfile);
}


/*
** set_ras_parms(): Set the RAS parameters or their conventional equivalents
**   (i.e. fermi level, etc).
**
*/
void set_ras_parms(void)
{
   int i,j,cnt;
   int errcod;
   int tot_expl_el,nras2alp,nras2bet,betsocc;
   int *ras1, *ras2, *ras3;
   int *orbsym;

   /* If the user asked for FCI=true, then override the other keywords
      if necessary to ensure that it's really a FCI
    */
   if (Parameters.fci == 1 && 
       (CalcInfo.num_alp_expl + CalcInfo.num_bet_expl) > Parameters.ex_lvl) 
   {
     Parameters.val_ex_lvl = 0;
     Parameters.ex_lvl = CalcInfo.num_alp_expl + CalcInfo.num_bet_expl;
     free(Parameters.ex_allow);
     Parameters.ex_allow = init_int_array(Parameters.ex_lvl);
     for (i=0; i<Parameters.ex_lvl; i++) Parameters.ex_allow[i] = 1;

     if (Parameters.print_lvl) {
       fprintf(outfile, "Note: Calculation requested is a full CI.\n");
       fprintf(outfile, 
               "Resetting EX_LVL to %d and turning on all excitations\n\n",
               Parameters.ex_lvl);
     }

   } /* end FCI override */

   /* reset ex_lvl if incompatible with number of electrons */
   if (Parameters.cc && (Parameters.cc_ex_lvl > 
     CalcInfo.num_alp_expl + CalcInfo.num_bet_expl)) {
     Parameters.cc_ex_lvl = CalcInfo.num_alp_expl + CalcInfo.num_bet_expl;
   }
   if (Parameters.ex_lvl > CalcInfo.num_alp_expl + CalcInfo.num_bet_expl) {
     Parameters.ex_lvl = CalcInfo.num_alp_expl + CalcInfo.num_bet_expl;
   }

   for (i=0,j=0; i<CalcInfo.nirreps; i++) j += CalcInfo.ras_opi[0][i];
   Parameters.a_ras1_lvl = Parameters.b_ras1_lvl = Parameters.ras1_lvl = j-1;

   /* figure out how many electrons are in RAS II */
   /* alpha electrons */
   for (i=0,nras2alp=0,betsocc=0; i<CalcInfo.nirreps; i++) {
      j = CalcInfo.docc[i] - CalcInfo.frozen_docc[i] - CalcInfo.ras_opi[0][i];
      if (Parameters.opentype == PARM_OPENTYPE_HIGHSPIN) {
         j += CalcInfo.socc[i];
         }
      else if (Parameters.opentype == PARM_OPENTYPE_SINGLET) {
         if (betsocc + CalcInfo.socc[i] <= CalcInfo.spab)
            betsocc += CalcInfo.socc[i];
         else {
            j += CalcInfo.socc[i] - (CalcInfo.spab - betsocc);
            betsocc = CalcInfo.spab;
            }
         }
      if (j > 0) nras2alp += j;
      if (j > CalcInfo.ras_opi[1][i]) {
         fprintf(outfile, "(set_ras_parms): detecting %d electrons ",
            j - CalcInfo.ras_opi[1][i]);
         fprintf(outfile, "in RAS III for irrep %d.\n", i);
         fprintf(outfile, "Some parts of DETCI assume all elec in I and II\n");
         }
      }
   /* beta electrons */
   for (i=0,nras2bet=0,betsocc=0; i<CalcInfo.nirreps; i++) {
      j = CalcInfo.docc[i] - CalcInfo.frozen_docc[i] - CalcInfo.ras_opi[0][i];
      if (Parameters.opentype == PARM_OPENTYPE_SINGLET && CalcInfo.socc[i]) {
         if (betsocc + CalcInfo.socc[i] <= CalcInfo.spab)
            j += CalcInfo.socc[i];
         else {
            j += CalcInfo.spab - betsocc;
            betsocc = CalcInfo.spab;
            }
         }
      if (j > 0) nras2bet += j;
      if (j > CalcInfo.ras_opi[1][i]) {
         fprintf(outfile, "(set_ras_parms): detecting %d electrons ",
            j - CalcInfo.ras_opi[1][i]);
         fprintf(outfile, "in RAS III for irrep %d.\n", i);
         fprintf(outfile, "Some parts of DETCI assume all elec in I and II\n");
         }
      }

   Parameters.a_ras1_max = (CalcInfo.num_alp_expl >
         Parameters.a_ras1_lvl + 1) ? Parameters.a_ras1_lvl + 1 :
         (CalcInfo.num_alp_expl) ;
   if (Parameters.fzc) Parameters.a_ras1_max += CalcInfo.num_fzc_orbs;

   Parameters.b_ras1_max = (CalcInfo.num_bet_expl >
         Parameters.b_ras1_lvl + 1) ? Parameters.b_ras1_lvl + 1:
         (CalcInfo.num_bet_expl) ;
   if (Parameters.fzc) Parameters.b_ras1_max += CalcInfo.num_fzc_orbs;

   for (i=0,j=0; i<CalcInfo.nirreps; i++) j += CalcInfo.ras_opi[1][i];
   Parameters.ras3_lvl = Parameters.ras1_lvl + j + 1;

   for (i=0,j=0; i<CalcInfo.nirreps; i++) j += CalcInfo.ras_opi[2][i];
   Parameters.ras4_lvl = Parameters.ras3_lvl + j;


   /* check Parameters to make sure everything consistent */

   if (Parameters.cc) {
     if (Parameters.cc_a_val_ex_lvl == -1)
       Parameters.cc_a_val_ex_lvl = Parameters.cc_val_ex_lvl;
     if (Parameters.cc_b_val_ex_lvl == -1)
       Parameters.cc_b_val_ex_lvl = Parameters.cc_val_ex_lvl;
     if (Parameters.cc_a_val_ex_lvl > Parameters.ras3_lvl - 
         Parameters.ras1_lvl - 1)
       Parameters.cc_a_val_ex_lvl = Parameters.ras3_lvl-Parameters.ras1_lvl-1;
     if (Parameters.cc_b_val_ex_lvl > Parameters.ras3_lvl - 
         Parameters.ras1_lvl - 1)
       Parameters.cc_b_val_ex_lvl = Parameters.ras3_lvl-Parameters.ras1_lvl-1;
     if (Parameters.cc_val_ex_lvl > Parameters.cc_a_val_ex_lvl +
         Parameters.cc_b_val_ex_lvl)
       Parameters.cc_val_ex_lvl = Parameters.cc_a_val_ex_lvl +
         Parameters.cc_b_val_ex_lvl;
   }

   /* deduce Parameters.cc_a_ras3_max and Parameters.cc_b_ras3_max if needed */
   if (Parameters.cc & (Parameters.cc_a_ras3_max == -1 || 
       Parameters.cc_b_ras3_max == -1)) {
      if (Parameters.cc_ras3_max != -1) { /* have parsed cc_ras3_max */
         Parameters.cc_a_ras3_max = 
            (Parameters.cc_ras3_max <= CalcInfo.num_alp_expl)
            ? Parameters.cc_ras3_max : CalcInfo.num_alp_expl;
         Parameters.cc_b_ras3_max = 
            (Parameters.cc_ras3_max <= CalcInfo.num_bet_expl)
            ? Parameters.cc_ras3_max : CalcInfo.num_bet_expl;
         }
      else {
         Parameters.cc_a_ras3_max = 
            (Parameters.cc_ex_lvl <= CalcInfo.num_alp_expl) 
            ? Parameters.cc_ex_lvl : CalcInfo.num_alp_expl;
         Parameters.cc_b_ras3_max = 
            (Parameters.cc_ex_lvl <= CalcInfo.num_bet_expl) 
            ? Parameters.cc_ex_lvl : CalcInfo.num_bet_expl; 
         }
      }

   if (Parameters.cc) {
      Parameters.a_ras3_max = 
         (Parameters.cc_a_ras3_max+2<=CalcInfo.num_alp_expl) 
         ? Parameters.cc_a_ras3_max+2 : CalcInfo.num_alp_expl;
      Parameters.b_ras3_max = 
         (Parameters.cc_b_ras3_max+2<=CalcInfo.num_bet_expl) 
         ? Parameters.cc_b_ras3_max+2 : CalcInfo.num_bet_expl;
      }

   if (Parameters.a_ras3_max == -1 || Parameters.b_ras3_max == -1) {
      if (Parameters.ras3_max != -1) { /* have parsed ras3_max */
         Parameters.a_ras3_max = (Parameters.ras3_max <= CalcInfo.num_alp_expl)
            ? Parameters.ras3_max : CalcInfo.num_alp_expl;
         Parameters.b_ras3_max = (Parameters.ras3_max <= CalcInfo.num_bet_expl)
            ? Parameters.ras3_max : CalcInfo.num_bet_expl;
         }
      else {
         Parameters.a_ras3_max = (Parameters.ex_lvl <= CalcInfo.num_alp_expl) 
            ? Parameters.ex_lvl : CalcInfo.num_alp_expl;
         Parameters.b_ras3_max = (Parameters.ex_lvl <= CalcInfo.num_bet_expl) 
            ? Parameters.ex_lvl : CalcInfo.num_bet_expl; 
         }
      }

   if (Parameters.cc) {
     if (Parameters.cc_ras4_max != -1) { /* have parsed */
        Parameters.cc_a_ras4_max = 
          (Parameters.cc_ras4_max <= CalcInfo.num_alp_expl)
          ? Parameters.cc_ras4_max : CalcInfo.num_alp_expl;
        Parameters.cc_b_ras4_max = 
          (Parameters.cc_ras4_max <= CalcInfo.num_bet_expl)
          ? Parameters.cc_ras4_max : CalcInfo.num_bet_expl;
        }
     else {
        Parameters.cc_a_ras4_max = Parameters.cc_a_ras3_max;
        Parameters.cc_b_ras4_max = Parameters.cc_b_ras3_max;
     }
     Parameters.a_ras4_max = 
       (Parameters.cc_a_ras4_max+2 <= CalcInfo.num_alp_expl) 
       ? Parameters.cc_a_ras4_max+2 : CalcInfo.num_alp_expl;
     Parameters.b_ras4_max = 
       (Parameters.cc_b_ras4_max+2 <= CalcInfo.num_bet_expl) 
       ? Parameters.cc_b_ras4_max+2 : CalcInfo.num_bet_expl;
   }
   else {
     if (Parameters.ras4_max != -1) { /* have parsed */
        Parameters.a_ras4_max = (Parameters.ras4_max <= CalcInfo.num_alp_expl)
          ? Parameters.ras4_max : CalcInfo.num_alp_expl;
        Parameters.b_ras4_max = (Parameters.ras4_max <= CalcInfo.num_bet_expl)
          ? Parameters.ras4_max : CalcInfo.num_bet_expl;
        }
     else {
        Parameters.a_ras4_max = Parameters.a_ras3_max;
        Parameters.b_ras4_max = Parameters.b_ras3_max;
     }
   }


   if (Parameters.cc) {
     if (Parameters.cc_ras34_max != -1) { /* have parsed */
       Parameters.cc_a_ras34_max = Parameters.cc_ras34_max;
       Parameters.cc_b_ras34_max = Parameters.cc_ras34_max;
     }
     else {
       Parameters.cc_a_ras34_max = Parameters.cc_a_ras3_max
                                 + Parameters.cc_a_ras4_max;
       Parameters.cc_b_ras34_max = Parameters.cc_b_ras3_max
                                 + Parameters.cc_b_ras4_max;
     }
     if (Parameters.ras34_max != -1) { /* have parsed */
        Parameters.a_ras34_max = Parameters.ras34_max;
        Parameters.b_ras34_max = Parameters.ras34_max;
        }
     else {
        Parameters.a_ras34_max = Parameters.cc_a_ras34_max+2;
        if (Parameters.a_ras34_max > CalcInfo.num_alp_expl)
          Parameters.a_ras34_max = CalcInfo.num_alp_expl;
        Parameters.b_ras34_max = Parameters.cc_b_ras34_max+2;
        if (Parameters.b_ras34_max > CalcInfo.num_bet_expl)
          Parameters.b_ras34_max = CalcInfo.num_bet_expl;
        Parameters.ras34_max = Parameters.cc_ras34_max+2;
        if (Parameters.ras34_max > Parameters.a_ras34_max +
            Parameters.b_ras34_max)
          Parameters.ras34_max = Parameters.a_ras34_max + 
            Parameters.b_ras34_max;
        }
   }
   else { /* non-CC */
     if (Parameters.ras34_max != -1) { /* have parsed */
       Parameters.a_ras34_max = Parameters.ras34_max;
       Parameters.b_ras34_max = Parameters.ras34_max;
     }
     else {
       Parameters.a_ras34_max = Parameters.a_ras3_max;
       Parameters.b_ras34_max = Parameters.b_ras3_max;
     }
   }

   i = Parameters.ras4_lvl - Parameters.ras3_lvl;
   if (Parameters.a_ras3_max > i) Parameters.a_ras3_max = i;
   if (Parameters.b_ras3_max > i) Parameters.b_ras3_max = i;
   if (Parameters.cc) {
     if (Parameters.cc_a_ras3_max > i) Parameters.cc_a_ras3_max = i;
     if (Parameters.cc_b_ras3_max > i) Parameters.cc_b_ras3_max = i;
   }

   i = CalcInfo.num_ci_orbs - Parameters.ras4_lvl;
   if (Parameters.a_ras4_max > i) Parameters.a_ras4_max = i;
   if (Parameters.b_ras4_max > i) Parameters.b_ras4_max = i;
   if (Parameters.cc) {
     if (Parameters.cc_a_ras4_max > i) Parameters.cc_a_ras4_max = i;
     if (Parameters.cc_b_ras4_max > i) Parameters.cc_b_ras4_max = i;
   }

   i = CalcInfo.num_ci_orbs - Parameters.ras3_lvl;
   if (Parameters.a_ras34_max > i) Parameters.a_ras34_max = i;
   if (Parameters.b_ras34_max > i) Parameters.b_ras34_max = i;
   if (Parameters.cc) {
     if (Parameters.cc_a_ras34_max > i) Parameters.cc_a_ras34_max = i;
     if (Parameters.cc_b_ras34_max > i) Parameters.cc_b_ras34_max = i;
   }

   i = (CalcInfo.num_alp_expl <= Parameters.a_ras1_lvl + 1) ? 
      CalcInfo.num_alp_expl : Parameters.a_ras1_lvl + 1;
   Parameters.a_ras1_min = i - Parameters.ex_lvl -
      Parameters.val_ex_lvl;
   if (Parameters.a_ras1_min < 0) Parameters.a_ras1_min = 0;
   Parameters.a_ras1_min += CalcInfo.num_fzc_orbs;
   Parameters.a_ras1_min += CalcInfo.num_cor_orbs;

   i = (CalcInfo.num_bet_expl <= Parameters.b_ras1_lvl + 1) ? 
      CalcInfo.num_bet_expl : Parameters.b_ras1_lvl + 1;
   Parameters.b_ras1_min = i - Parameters.ex_lvl -
      Parameters.val_ex_lvl;
   if (Parameters.b_ras1_min < 0) Parameters.b_ras1_min = 0;
   Parameters.b_ras1_min += CalcInfo.num_fzc_orbs;
   Parameters.b_ras1_min += CalcInfo.num_cor_orbs;

   tot_expl_el = CalcInfo.num_alp_expl + CalcInfo.num_bet_expl;
   if (Parameters.cc) {
     if (Parameters.cc_val_ex_lvl != 0) i = Parameters.cc_val_ex_lvl;
     else i = Parameters.cc_ex_lvl;
     if (Parameters.cc_ras3_max == -1) {
       Parameters.cc_ras3_max = (i <= tot_expl_el) ?  i : tot_expl_el ;
     }
     else {
       if (Parameters.cc_ras3_max > tot_expl_el) 
         Parameters.cc_ras3_max = tot_expl_el;
     } 
     if (Parameters.ras3_max == -1) 
       Parameters.ras3_max = Parameters.cc_ras3_max + 2;
   }
   if (Parameters.ras3_max == -1) {
      Parameters.ras3_max = (Parameters.ex_lvl <= tot_expl_el) ?
         Parameters.ex_lvl : tot_expl_el ;
      }
   else {
      if (Parameters.ras3_max > tot_expl_el) 
         Parameters.ras3_max = tot_expl_el;
      }

   i = 2 * (Parameters.ras4_lvl - Parameters.ras3_lvl);
   if (i < Parameters.ras3_max) Parameters.ras3_max = i;
   if (Parameters.cc) {
     if (i < Parameters.cc_ras3_max) Parameters.cc_ras3_max = i;
   }


   i = (tot_expl_el < 2*(Parameters.ras1_lvl + 1)) ? tot_expl_el :
      2*(Parameters.ras1_lvl + 1) ;
   
   Parameters.ras1_min = i - Parameters.ex_lvl - 
      Parameters.val_ex_lvl + 2 * CalcInfo.num_fzc_orbs;
   
   if (Parameters.a_ras1_min + Parameters.b_ras1_min > Parameters.ras1_min)
      Parameters.ras1_min = Parameters.a_ras1_min + Parameters.b_ras1_min;

   if (Parameters.cc && Parameters.cc_ras4_max == -1) {
     Parameters.cc_ras4_max = (Parameters.cc_ex_lvl <= tot_expl_el) ?
       Parameters.cc_ex_lvl : tot_expl_el;
   }

   if (Parameters.ras4_max == -1) {
     if (Parameters.cc) {
       Parameters.ras4_max = (Parameters.cc_ras4_max+2 <= tot_expl_el) ?
         Parameters.cc_ras4_max+2 : tot_expl_el;
     }
     else Parameters.ras4_max = (Parameters.ex_lvl <= tot_expl_el) ?
       Parameters.ex_lvl : tot_expl_el;
   }

   i = 2 * (CalcInfo.num_ci_orbs - Parameters.ras4_lvl);
   if (i < Parameters.ras4_max) Parameters.ras4_max = i;
   if (Parameters.cc) {
     if (i < Parameters.cc_ras4_max) Parameters.cc_ras4_max = i;
   }

   if (Parameters.cc && Parameters.cc_ras34_max == -1)
     Parameters.cc_ras34_max = Parameters.cc_ras3_max + Parameters.cc_ras4_max;
   i = 2 * (CalcInfo.num_ci_orbs - Parameters.ras3_lvl);
   if (Parameters.cc) {
     if (i < Parameters.cc_ras34_max) Parameters.cc_ras34_max = i;
   }

   if (Parameters.ras34_max == -1 && !Parameters.cc) 
     Parameters.ras34_max = Parameters.ras3_max;
   else
     Parameters.ras34_max = Parameters.cc_ras34_max + 2;

   i = 2 * (CalcInfo.num_ci_orbs - Parameters.ras3_lvl);
   if (i < Parameters.ras34_max) Parameters.ras34_max = i;
   if (Parameters.ras34_max > tot_expl_el) Parameters.ras34_max = tot_expl_el;

   if (Parameters.a_ras34_max > Parameters.a_ras3_max + Parameters.a_ras4_max)
     Parameters.a_ras34_max = Parameters.a_ras3_max + Parameters.a_ras4_max;
   if (Parameters.cc && (Parameters.cc_a_ras34_max > 
     Parameters.cc_a_ras3_max + Parameters.cc_b_ras4_max))
     Parameters.cc_a_ras34_max = Parameters.cc_a_ras3_max + 
     Parameters.cc_a_ras4_max;

   if (Parameters.b_ras34_max > Parameters.b_ras3_max + Parameters.b_ras4_max)
     Parameters.b_ras34_max = Parameters.b_ras3_max + Parameters.b_ras4_max;
   if (Parameters.cc && (Parameters.cc_b_ras34_max > 
     Parameters.cc_b_ras3_max + Parameters.cc_b_ras4_max))
     Parameters.cc_b_ras34_max = Parameters.cc_b_ras3_max + 
     Parameters.cc_b_ras4_max;

   /* now just re-check some basic things */
   if (Parameters.a_ras34_max > CalcInfo.num_alp_expl)
     Parameters.a_ras34_max = CalcInfo.num_alp_expl;
   if (Parameters.b_ras34_max > CalcInfo.num_bet_expl)
     Parameters.b_ras34_max = CalcInfo.num_bet_expl;
   if (Parameters.cc) {
     if (Parameters.cc_a_ras34_max > CalcInfo.num_alp_expl)
       Parameters.cc_a_ras34_max = CalcInfo.num_alp_expl;
     if (Parameters.cc_b_ras34_max > CalcInfo.num_bet_expl)
       Parameters.cc_b_ras34_max = CalcInfo.num_bet_expl;
   }

   if (Parameters.ras34_max > Parameters.a_ras34_max + Parameters.b_ras34_max)
     Parameters.ras34_max = Parameters.a_ras34_max + Parameters.b_ras34_max;

}


/*
** print_ras_parms(): Set the RAS parameters or their conventional equivalents
**   (i.e. fermi level, etc).
**
*/
void print_ras_parms(void)
{
  int i, j;

  fprintf(outfile, "ORBITALS:\n") ;
  fprintf(outfile, "   NMO          =   %6d      NUM ALP      =   %6d\n",
    CalcInfo.nmo, CalcInfo.num_alp);
  fprintf(outfile, "   ORBS IN CI   =   %6d      NUM ALP EXPL =   %6d\n",
    CalcInfo.num_ci_orbs, CalcInfo.num_alp_expl);
  fprintf(outfile, "   FROZEN CORE  =   %6d      NUM BET      =   %6d\n",
    CalcInfo.num_fzc_orbs, CalcInfo.num_bet);
  fprintf(outfile, "   RESTR CORE   =   %6d      NUM BET EXPL =   %6d\n",
    CalcInfo.num_cor_orbs, CalcInfo.num_bet_expl);
  fprintf(outfile, "   IOPEN        =   %6s\n", CalcInfo.iopen ? "yes" :
    "no");
  fprintf(outfile, "   RAS1 LVL     =   %6d      A RAS3 MAX   =   %6d\n",
    Parameters.ras1_lvl, Parameters.a_ras3_max);
  fprintf(outfile, "   RAS1 MIN     =   %6d      B RAS3 MAX   =   %6d\n",
    Parameters.ras1_min, Parameters.b_ras3_max);
  fprintf(outfile, "   A RAS1 LVL   =   %6d      RAS4 LVL     =   %6d\n", 
    Parameters.a_ras1_lvl, Parameters.ras4_lvl);
  fprintf(outfile, "   A RAS1 MIN   =   %6d      A RAS4 MAX   =   %6d\n", 
    Parameters.a_ras1_min, Parameters.a_ras4_max);
  fprintf(outfile, "   A RAS1 MAX   =   %6d      B RAS4 MAX   =   %6d\n", 
    Parameters.a_ras1_max, Parameters.b_ras4_max);
  fprintf(outfile, "   B RAS1 LVL   =   %6d      RAS4 MAX     =   %6d\n", 
    Parameters.b_ras1_lvl, Parameters.ras4_max);
  fprintf(outfile, "   B RAS1 MIN   =   %6d      A RAS34 MAX  =   %6d\n", 
    Parameters.b_ras1_min, Parameters.a_ras34_max);
  fprintf(outfile, "   B RAS1 MAX   =   %6d      B RAS34 MAX  =   %6d\n", 
    Parameters.b_ras1_max, Parameters.b_ras34_max);
  fprintf(outfile, "   RAS3 LVL     =   %6d      RAS34 MAX    =   %6d\n", 
    Parameters.ras3_lvl, Parameters.ras34_max);
  fprintf(outfile, "   RAS3 MAX     =   %6d\n", Parameters.ras3_max);
  if (Parameters.cc) {
    fprintf(outfile, "   CC RAS3 MAX  =   %6d      CC RAS4 MAX  =   %6d\n",
      Parameters.cc_ras3_max, Parameters.cc_ras4_max);
    fprintf(outfile, "   CC A RAS3 MAX=   %6d      CC B RAS3 MAX=   %6d\n",
      Parameters.cc_a_ras3_max, Parameters.cc_b_ras3_max);
    fprintf(outfile, "   CC A RAS4 MAX=   %6d      CC B RAS4 MAX=   %6d\n", 
      Parameters.cc_a_ras4_max, Parameters.cc_b_ras4_max);
    fprintf(outfile, "   CC RAS34 MAX =   %6d\n",
      Parameters.cc_ras34_max);
    fprintf(outfile, "   CC A RAS34 MAX=  %6d      CC B RAS34 MAX=  %6d\n", 
      Parameters.cc_a_ras34_max, Parameters.cc_b_ras34_max);
    fprintf(outfile, "   CC MIXED     =   %6s      CC FIX EXTERN =  %6s\n",
      Parameters.cc_mixed ? "yes" : "no",
      Parameters.cc_fix_external ? "yes" : "no");
    fprintf(outfile, "   CC VARIATIONAL=  %6s\n",
      Parameters.cc_variational ? "yes" : "no");
  }

  fprintf(outfile, "\n");
  fprintf(outfile, "   DOCC         = ") ;
  for (i=0; i<CalcInfo.nirreps; i++) {
    fprintf(outfile, "%2d ", CalcInfo.docc[i]) ;
  }
  fprintf(outfile, "\n   SOCC         = ") ;
  for (i=0; i<CalcInfo.nirreps; i++) {
    fprintf(outfile, "%2d ", CalcInfo.socc[i]) ;
  }
  fprintf(outfile, "\n   FROZEN DOCC  = ") ;
  for (i=0; i<CalcInfo.nirreps; i++) {
    fprintf(outfile, "%2d ", CalcInfo.frozen_docc[i]) ;
  }
  fprintf(outfile, "\n   FROZEN UOCC  = ") ;
  for (i=0; i<CalcInfo.nirreps; i++) {
    fprintf(outfile, "%2d ", CalcInfo.frozen_uocc[i]) ;
  }
  fprintf(outfile, "\n");
  for (i=0; i<4; i++) {
    fprintf(outfile, "   RAS %d        = ",i+1);
    for (j=0; j<CalcInfo.nirreps; j++) {
      fprintf(outfile,"%2d ",CalcInfo.ras_opi[i][j]);
    }
    fprintf(outfile, "\n");
  }

  fprintf(outfile,
     "*******************************************************\n\n");
}

}} // namespace psi::detci

