/*! \file
    \ingroup OPTKING
    \brief GET_OPTINFO   reads optimization parameters from input.dat
*/

#define EXTERN
#include "globals.h"
#undef EXTERN
#include "cartesians.h"
#include "simples.h"
#include "salc.h"
#include "opt.h"

#include <libqt/qt.h>
#include <libipv1/ip_lib.h>
#include <libpsio/psio.h>
#include <libchkpt/chkpt.h>

namespace psi { //namespace optking {

static double power(double x, int y);

// opt_mode == name of calling function/ name of function called from psi4
void get_optinfo(void) {

  static bool read_options = false;

  if (read_options) { // already done
    intro("options read previously");  
    return;
  }
  
  read_options = true;
  intro();

  ip_cwk_add(":OPTKING");

  // determine if simples and salcs are provided in intco.dat 
  optinfo.simples_present = optinfo.salcs_present = false;
  opt_ffile_noexit(&fp_intco, "intco.dat", 2);
  if (fp_intco != NULL)
    ip_append(fp_intco, outfile) ;
  if ( ip_exist(":INTCO",0) ) {
    ip_cwk_add(":INTCO");
    optinfo.simples_present = true;
  }
  if ( ip_exist("SYMM",0) || ip_exist("ASYMM",0) )
    optinfo.salcs_present = true;
  if (fp_intco != NULL)
    fclose(fp_intco);

  optinfo.iteration = 0;
  optinfo.micro_iteration = 0;

  open_PSIF();
  if (psio_tocscan(PSIF_OPTKING, "Iteration") != NULL)
    psio_read_entry(PSIF_OPTKING, "Iteration", (char *) &(optinfo.iteration),sizeof(int));
  if (psio_tocscan(PSIF_OPTKING, "Micro_iteration") != NULL)
    psio_read_entry(PSIF_OPTKING, "Micro_iteration", (char *) &(optinfo.micro_iteration),sizeof(int));
  if (psio_tocscan(PSIF_OPTKING, "Balked last time") != NULL)
    psio_read_entry(PSIF_OPTKING, "Balked last time", (char *) &(optinfo.balked_last_time), sizeof(int));
  close_PSIF();

  char *junk;
  int a, i, cnt, cnt2, natom, nallatom, errcod;
  double tval;

  optinfo.dertype = 0;
  if (ip_exist("DERTYPE",0)) {
    errcod = ip_string("DERTYPE", &(junk),0);
    if (errcod != IPE_OK) optinfo.dertype = 0;
    else if(!strcmp(junk,"NONE")) optinfo.dertype = 0;
    else if(!strcmp(junk,"FIRST")) optinfo.dertype = 1;
    else if(!strcmp(junk,"SECOND")) optinfo.dertype = 2;
    else if(!strcmp(junk,"RESPONSE")) optinfo.dertype = 3; /* linear response */
    else {
      printf("Invalid value of input keyword DERTYPE: %s\n", junk);
      exit(PSI_RETURN_FAILURE);
    }
    free(junk);
  }

  errcod = ip_string("WFN", &(optinfo.wfn), 0);
  errcod = ip_string("JOBTYPE", &(optinfo.jobtype), 0);
 
  optinfo.numerical_dertype = 0;
  if (ip_exist("DERTYPE",0)) {
    if ( (strcmp(optinfo.jobtype,"OPT")==0) && (optinfo.dertype==0) ) 
      optinfo.numerical_dertype = 1;
  }

  // only use listed force constants for finite differences
  if (ip_exist("SELECTED_FC",0))
    optinfo.selected_fc = true;

  optinfo.nonselected_fc = OPTInfo::EMPIRICAL;
  if (ip_exist("NONSELECTED_FC",0)) {
    errcod = ip_string("NONSELECTED_FC",&(junk),0);
    if (!strcmp(junk,"KEEP"))
      optinfo.nonselected_fc = OPTInfo::KEEP;
    else if (!strcmp(junk,"EMPIRICAL"))
      optinfo.nonselected_fc = OPTInfo::EMPIRICAL;
    else
      fprintf(outfile,"Unable to understand NONSELECTED_FC keyword entry.\n");
    free(junk);
  }

  /* optimize in cartesian coordinates and ignore redundant/delocalize keywords */
  optinfo.cartesian = 0;
  ip_boolean("CARTESIAN", &(optinfo.cartesian),0);

  optinfo.redundant = 1; optinfo.delocalize = 0;

  if ( optinfo.numerical_dertype == 1 )
    { optinfo.redundant = 0; optinfo.delocalize =1; }

  /// if user selects fc, assume the ones present are redundant or OK anyway
  if ((optinfo.mode == MODE_DISP_IRREP) || (optinfo.mode == MODE_DISP_NOSYMM) ) {
    if (!optinfo.selected_fc) {
      optinfo.redundant = 0; optinfo.delocalize =1;
    }
  }

  ip_boolean("DELOCALIZE", &(optinfo.delocalize),0);
  if (optinfo.delocalize)
    optinfo.redundant = 0;

  /* default for fragments is to use R */
  optinfo.frag_dist_rho = 0;
  ip_boolean("FRAGMENT_DISTANCE_INVERSE", &(optinfo.frag_dist_rho),0);
  optinfo.fix_intrafragment = 0;
  ip_boolean("FIX_INTRAFRAGMENT", &(optinfo.fix_intrafragment),0);
  optinfo.fix_interfragment = 0;
  ip_boolean("FIX_INTERFRAGMENT", &(optinfo.fix_interfragment),0);
  optinfo.analytic_interfragment = 1;
  ip_boolean("ANALYTIC_INTERFRAGMENT", &(optinfo.analytic_interfragment),0);

  /* print options */
  optinfo.print_simples = 0;
  ip_boolean("PRINT_SIMPLES", &(optinfo.print_simples),0);
  optinfo.print_params = 0;
  ip_boolean("PRINT_PARAMS", &(optinfo.print_params),0);
  optinfo.print_delocalize = 0;
  ip_boolean("PRINT_DELOCALIZE", &(optinfo.print_delocalize),0);
  optinfo.print_symmetry = 0;
  ip_boolean("PRINT_SYMMETRY", &(optinfo.print_symmetry),0);
  optinfo.print_hessian = 0;
  ip_boolean("PRINT_HESSIAN", &(optinfo.print_hessian),0);
  optinfo.print_cartesians = 0;
  ip_boolean("PRINT_CARTESIANS", &(optinfo.print_cartesians),0);
  optinfo.print_fconst = 0;
  ip_boolean("PRINT_FCONST", &(optinfo.print_fconst),0);
  optinfo.print_debug_backtransformation = 0;
  ip_boolean("PRINT_DEBUG_BACKTRANSFORMATION", &(optinfo.print_debug_backtransformation),0);

  /* optimization parameters */
  optinfo.optimize = 1;
  if (ip_exist("DISPLACEMENTS",0)) {
    if ((optinfo.mode != MODE_DISP_LOAD) && (optinfo.mode != MODE_ENERGY_SAVE)
      && (optinfo.mode != MODE_GRAD_ENERGY)) {
      optinfo.do_only_displacements = 1;
      optinfo.mode = MODE_DISP_USER;
      optinfo.optimize = 0;
    }
  }

  optinfo.freq_irrep = -1;
  if (ip_exist("FREQ_IRREP",0)) {
    ip_data("FREQ_IRREP","%d",&(optinfo.freq_irrep),0);
  }

  optinfo.H_update = OPTInfo::BFGS;
  if (ip_exist("H_UPDATE",0)) {
    errcod = ip_string("H_UPDATE",&(junk),0);
    if (!strcmp(junk,"NONE"))
      optinfo.H_update = OPTInfo::NONE;
    else if (!strcmp(junk,"BFGS"))
      optinfo.H_update = OPTInfo::BFGS;
    else if (!strcmp(junk,"MS"))
      optinfo.H_update = OPTInfo::MS;
    else if (!strcmp(junk,"POWELL"))
      optinfo.H_update = OPTInfo::POWELL;
    else if (!strcmp(junk,"BOFILL"))
      optinfo.H_update = OPTInfo::BOFILL;
    else
      fprintf(outfile,"Unable to understand H_UPDATE keyword entry.\n");
    free(junk);
  }

  optinfo.empirical_H = OPTInfo::SCHLEGEL;
  if (ip_exist("EMPIRICAL_H",0)) {
    errcod = ip_string("EMPIRICAL_H",&(junk),0);
    if (!strcmp(junk,"FISCHER"))
      optinfo.empirical_H = OPTInfo::FISCHER;
    else if (!strcmp(junk,"SCHLEGEL"))
      optinfo.empirical_H = OPTInfo::SCHLEGEL;
    else
      fprintf(outfile,"Unable to understand EMPIRICAL_H keyword entry.\n");
    free(junk);
  }

  i=0;
  ip_boolean("RFO", &i,0);
  optinfo.rfo = (i ? true : false);

  i=0;
  ip_boolean("RFO_FOLLOW_ROOT", &i,0);
  optinfo.rfo_follow_root = (i ? true : false);

  i=1;
  ip_data("RFO_ROOT", "%d", &i,0);
  optinfo.rfo_root = i-1;

  i=3;
  ip_data("MAX_CONSECUTIVE_LINE_SEARCHES", "%d", &i,0);
  optinfo.max_consecutive_line_searches = i;

  i = 8;
  ip_data("LINE_SEARCH_MIN","%d",&i,0);
  optinfo.line_search_min = power(10.0,-1*i);

  optinfo.step_energy_limit = 0.4; // fraction error in energy prediction to tolerate
  errcod = ip_data("STEP_ENERGY_LIMIT", "%lf", &tval,0);
  if (errcod == IPE_OK) optinfo.step_energy_limit = tval;

  //fraction error to accept in guess step backward
  optinfo.step_energy_limit_back = 2.0/3.0*optinfo.step_energy_limit;
  errcod = ip_data("STEP_ENERGY_LIMIT_BACK", "%lf", &tval,0);
  if (errcod == IPE_OK) optinfo.step_energy_limit_back = tval;

  // ACS (11/06) Are we getting our energies from outside of PSI3?
  optinfo.external_energies = 0;
  ip_boolean("EXTERNAL_ENERGIES",&(optinfo.external_energies),0);

  optinfo.H_update_use_last = 3;
  ip_data("H_UPDATE_USE_LAST","%d",&(optinfo.H_update_use_last),0);
  optinfo.mix_types = 1;
  ip_boolean("MIX_TYPES",&(optinfo.mix_types),0);

  ip_data("POINTS","%d",&(optinfo.points),0);
  /*
  if ( (optinfo.dertype == 0) && (optinfo.numerical_dertype == 2) ) {
    fprintf(outfile,"O(h^2) formula used for frequencies by energy points.\n");
    optinfo.points = 3;
  }
  */

  optinfo.step_limit_cart = 0.3;
  errcod = ip_data("STEP_LIMIT_CART","%lf",&tval,0);
  if (errcod == IPE_OK)
    optinfo.step_limit_cart = tval;

  // max change in a simple internal coordinate in bohr or radian
  optinfo.step_limit = 0.4;
  errcod = ip_data("STEP_LIMIT","%lf",&tval,0);
  if (errcod == IPE_OK)
    optinfo.step_limit = tval;

  /* takes values of 1,2,3 for x,y,z for location of first dummy of linear bend*/
  optinfo.dummy_axis_1 = 1;
  ip_data("DUMMY_AXIS_1","%d",&(optinfo.dummy_axis_1),0);
  optinfo.dummy_axis_1 -= 1;

  optinfo.zmat = 0;
  if (ip_exist("ZMAT",0)) optinfo.zmat = 1;

  optinfo.zmat_simples = 0;
  ip_boolean("ZMAT_SIMPLES",&(optinfo.zmat_simples),0);

  optinfo.test_B = 0;
  ip_boolean("TEST_B",&(optinfo.test_B),0);

  if (ip_exist("OPT_CONV",0)) {
    errcod = ip_string("OPT_CONV", &junk, 0);
    if (errcod != IPE_OK) {
      fprintf(outfile,"OPT_CONV should be one of {LOOSE, NORMAL, TIGHT, VERY_TIGHT, BAKER, QCHEM, G03_NORMAL, G03_TIGHT, G03_VERY_TIGHT, GENERAL}\n");
      throw("Could not read OPT_CONV from input.");
    }
    else if (!strcmp(junk,"LOOSE")) {
      optinfo.opt_conv       = OPTInfo::LOOSE;
      optinfo.conv_max_force = power(10.0, -3);
      optinfo.conv_max_DE    = power(10.0, -6);
      optinfo.conv_max_disp  = power(10.0, -2);
    }
    else if (!strcmp(junk,"NORMAL")) { // also below
      optinfo.opt_conv       = OPTInfo::NORMAL;
      optinfo.conv_max_force = power(10.0, -4);
      optinfo.conv_max_DE    = power(10.0, -8);
      optinfo.conv_max_disp  = power(10.0, -3);
    }
    else if (!strcmp(junk,"TIGHT")) {
      optinfo.opt_conv       = OPTInfo::TIGHT;
      optinfo.conv_max_force = power(10.0, -5);
      optinfo.conv_max_DE    = power(10.0, -10);
      optinfo.conv_max_disp  = power(10.0, -4);
    }
    else if (!strcmp(junk,"VERY_TIGHT")) {
      optinfo.opt_conv       = OPTInfo::VERY_TIGHT;
      optinfo.conv_max_force = power(10.0, -6);
      optinfo.conv_max_DE    = power(10.0, -12);
      optinfo.conv_max_disp  = power(10.0, -5);
    }
    else if (!strcmp(junk,"BAKER")) {
      optinfo.opt_conv       = OPTInfo::BAKER;
      optinfo.conv_max_force = 3.0e-4;
      optinfo.conv_max_DE    = power(10.0, -6);
      optinfo.conv_max_disp  = 3.0e-4;
    }
    else if (!strcmp(junk,"QCHEM")) {
      optinfo.opt_conv       = OPTInfo::QCHEM;
      optinfo.conv_max_force = 3.0e-4;
      optinfo.conv_max_DE    = power(10.0, -6);
      optinfo.conv_max_disp  = 1.2e-3;
    }
    else if (!strcmp(junk,"G03_NORMAL")) {
      optinfo.opt_conv       = OPTInfo::G03_NORMAL;
      optinfo.conv_max_force = 4.5e-4;
      optinfo.conv_rms_force = 3.0e-4;
      optinfo.conv_max_disp  = 1.8e-3;
      optinfo.conv_rms_disp  = 1.2e-3;
    }
    else if (!strcmp(junk,"G03_TIGHT")) {
      optinfo.opt_conv       = OPTInfo::G03_TIGHT;
      optinfo.conv_max_force = 1.5e-5;
      optinfo.conv_rms_force = 1.0e-5;
      optinfo.conv_max_disp  = 6.0e-5;
      optinfo.conv_rms_disp  = 4.0e-5;
    }
    else if (!strcmp(junk,"G03_VERY_TIGHT")) {
      optinfo.opt_conv       = OPTInfo::G03_VERY_TIGHT;
      optinfo.conv_max_force = 2.0e-6;
      optinfo.conv_rms_force = 1.0e-6;
      optinfo.conv_max_disp  = 6.0e-6;
      optinfo.conv_rms_disp  = 4.0e-6;
    }
    else if (!strcmp(junk,"GENERAL")) {
      optinfo.opt_conv       = OPTInfo::GENERAL;
      optinfo.conv_max_force = 0.0;
      optinfo.conv_rms_force = 0.0;
      optinfo.conv_max_DE    = 0.0;
      optinfo.conv_max_disp  = 0.0;
      optinfo.conv_rms_disp  = 0.0;
    }
    else {
      fprintf(outfile,"OPT_CONV should be one of {LOOSE, NORMAL, TIGHT, VERY_TIGHT, BAKER, QCHEM, G03_NORMAL, G03_TIGHT, G03_VERY_TIGHT, GENERAL}\n");
      throw("Could not read OPT_CONV from input.");
    }
    free(junk);
  }
  else { // default NORMAL criteria
    optinfo.opt_conv       = OPTInfo::NORMAL;
    optinfo.conv_max_force = power(10.0, -4);
    optinfo.conv_max_DE    = power(10.0, -8);
    optinfo.conv_max_disp  = power(10.0, -3);
  }

  if ( (optinfo.opt_conv == OPTInfo::LOOSE) || (optinfo.opt_conv == OPTInfo::NORMAL) || (optinfo.opt_conv == OPTInfo::TIGHT) ||
       (optinfo.opt_conv == OPTInfo::VERY_TIGHT) || (optinfo.opt_conv == OPTInfo::BAKER) || (optinfo.opt_conv == OPTInfo::QCHEM) ) {

    if (ip_exist("CONV_MAX_FORCE",0)) {
      errcod = ip_data("CONV_MAX_FORCE","%lf",&tval,0);
      if (errcod == IPE_OK) optinfo.conv_max_force = tval;
    }

    if (ip_exist("CONV_MAX_DE",0)) {
      errcod = ip_data("CONV_MAX_DE","%lf",&tval,0);
      if (errcod == IPE_OK) optinfo.conv_max_DE = tval;
    }

    if (ip_exist("CONV_MAX_DISP",0)) {
      errcod = ip_data("CONV_MAX_DISP","%lf",&tval,0);
      if (errcod == IPE_OK) optinfo.conv_max_disp = tval;
    }
  }

  else if ( (optinfo.opt_conv == OPTInfo::G03_NORMAL) || (optinfo.opt_conv == OPTInfo::G03_TIGHT) || (optinfo.opt_conv == OPTInfo::G03_VERY_TIGHT) ) {

    if (ip_exist("CONV_MAX_FORCE",0)) {
      errcod = ip_data("CONV_MAX_FORCE","%lf",&tval,0);
      if (errcod == IPE_OK) optinfo.conv_max_force = tval;
    }

    if (ip_exist("CONV_RMS_FORCE",0)) {
      errcod = ip_data("CONV_RMS_FORCE","%lf",&tval,0);
      if (errcod == IPE_OK) optinfo.conv_rms_force = tval;
    }

    if (ip_exist("CONV_MAX_DISP",0)) {
      errcod = ip_data("CONV_MAX_DISP","%lf",&tval,0);
      if (errcod == IPE_OK) optinfo.conv_max_disp = tval;
    }

    if (ip_exist("CONV_RMS_DISP",0)) {
      errcod = ip_data("CONV_RMS_DISP","%lf",&tval,0);
      if (errcod == IPE_OK) optinfo.conv_rms_disp = tval;
    }
  }

  else if ( optinfo.opt_conv == OPTInfo::GENERAL ) {
    a = 0;

    if (ip_exist("CONV_MAX_FORCE",0)) {
      errcod = ip_data("CONV_MAX_FORCE","%lf",&tval,0);
      if (errcod == IPE_OK) optinfo.conv_max_force = tval;
      a++;
    }

    if (ip_exist("CONV_RMS_FORCE",0)) {
      errcod = ip_data("CONV_RMS_FORCE","%lf",&tval,0);
      if (errcod == IPE_OK) optinfo.conv_rms_force = tval;
      a++;
    }

    if (ip_exist("CONV_MAX_DE",0)) {
      errcod = ip_data("CONV_MAX_DE","%lf",&tval,0);
      if (errcod == IPE_OK) optinfo.conv_max_DE = tval;
      a++;
    }

    if (ip_exist("CONV_MAX_DISP",0)) {
      errcod = ip_data("CONV_MAX_DISP","%lf",&tval,0);
      if (errcod == IPE_OK) optinfo.conv_max_disp = tval;
      a++;
    }

    if (ip_exist("CONV_RMS_DISP",0)) {
      errcod = ip_data("CONV_RMS_DISP","%lf",&tval,0);
      if (errcod == IPE_OK) optinfo.conv_rms_disp = tval;
      a++;
    }

    if (a == 0) {
      fprintf(outfile,"OPT_CONV=GENERAL should be accompanied by at least one of {CONV_MAX_FORCE, CONV_RMS_FORCE, CONV_MAX_DE, CONV_MAX_DISP, CONV_RMS_DISP}\n");
      throw("Could not find convergence parameters for OPT_CONV=GENERAL in input.");
    }
  }

  a= 5;
  ip_data("EV_TOL","%d",&a,0);
  optinfo.ev_tol = power(10.0, -1*a);

  optinfo.scale_connectivity = 1.3;
  ip_data("SCALE_CONNECTIVITY","%lf",&(optinfo.scale_connectivity),0);

  optinfo.disp_size = 0.005;
  ip_data("EDISP","%lf",&(optinfo.disp_size),0);

  /* back-transformation parameters */
  optinfo.bt_max_iter = 25;
  ip_data("BT_MAX_ITER","%d",&(optinfo.bt_max_iter),0);
  a = 4; /* some minimal level of this seems necessary for butane */
  ip_data("BT_DQ_CONV","%d",&a,0);
  optinfo.bt_dq_conv  = power(10.0, -1*a);
  a = 8;
  ip_data("BT_DX_CONV","%d",&a,0);
  optinfo.bt_dx_conv  = power(10.0, -1*a);


  /* These determing how seriously to take torsions very near 180 */
  optinfo.cos_tors_near_1_tol    =  0.9999999999;
  ip_data("COS_TORS_NEAR_1_TOL","%lf",&(optinfo.cos_tors_near_1_tol),0);
  optinfo.cos_tors_near_neg1_tol = -0.9999999999;
  ip_data("COS_TORS_NEAR_NEG1_TOL","%lf",&(optinfo.cos_tors_near_neg1_tol),0);
  optinfo.sin_phi_denominator_tol = 1.0e-7;
  ip_data("SIN_PHI_DENOMINATOR_TOL","%lf",&(optinfo.sin_phi_denominator_tol),0);


  /* Compute dummy atom lookup arrays */

  chkpt_init(PSIO_OPEN_OLD);
  optinfo.natom = natom = chkpt_rd_natom();
  optinfo.nallatom = nallatom = chkpt_rd_nallatom();
  optinfo.atom_dummy = chkpt_rd_atom_dummy();
  optinfo.nfragment = chkpt_rd_nfragment();
  if (optinfo.nfragment > 1) {
    optinfo.natom_per_fragment = chkpt_rd_natom_per_fragment();
    optinfo.nallatom_per_fragment = chkpt_rd_nallatom_per_fragment();
    optinfo.nref_per_fragment = chkpt_rd_nref_per_fragment();
    optinfo.fragment_coeff = chkpt_rd_fragment_coeff();
  }
  chkpt_close();

  optinfo.to_dummy = init_int_array(natom);
  optinfo.to_nodummy = init_int_array(nallatom);

  cnt=0; cnt2=0;
  for (i=0; i<nallatom; ++i) {
    if (optinfo.atom_dummy[i]) continue;
    else optinfo.to_dummy[cnt++] = i;
    
    if (optinfo.atom_dummy[i]) continue;
    else optinfo.to_nodummy[i] = cnt2++;
  }

  if (optinfo.print_params) {
    for (i=0;i<nallatom;++i)
      fprintf(outfile,"atom_dummy[%d]: %d\n",i,optinfo.atom_dummy[i]);
    for (i=0;i<nallatom;++i)
      fprintf(outfile,"to_nodummy[%d]: %d\n",i,optinfo.to_nodummy[i]);
    for (i=0;i<natom;++i)
      fprintf(outfile,"to_dummy[%d]: %d\n",i,optinfo.to_dummy[i]);
    fflush(outfile);
  }

  /*
  fprintf(outfile,"\nIteration: %d     ",optinfo.iteration+1);
  if (optinfo.numerical_dertype > 0)
    fprintf(outfile,"Micro_iteration: %d",optinfo.micro_iteration+1);
  fprintf(outfile,"\n");
  */
  if (optinfo.print_params) {
    fprintf(outfile,"\n+++ Optinfo Parameters +++\n");
    fprintf(outfile,"print_params:  %d\n",optinfo.print_params);
    fprintf(outfile,"print_simples: %d\n",optinfo.print_simples);
    fprintf(outfile,"print_delocalize %d\n",optinfo.print_delocalize);
    fprintf(outfile,"print_symmetry %d\n",optinfo.print_symmetry);
    fprintf(outfile,"optimize:      %d\n",optinfo.optimize);
    fprintf(outfile,"zmat:          %d\n",optinfo.zmat);
    fprintf(outfile,"dummy_axis_1:    %d\n",optinfo.dummy_axis_1);
    fprintf(outfile,"sacc fd points:  %d\n",optinfo.points);
    fprintf(outfile,"zmat_simples:  %d\n",optinfo.zmat_simples);
    fprintf(outfile,"redundant:     %d\n",optinfo.redundant);
    fprintf(outfile,"cartesian:     %d\n",optinfo.cartesian);
    fprintf(outfile,"searching for: %d\n", optinfo.ts ? "transition state" : "minimum");

    if (optinfo.H_update == OPTInfo::NONE)
      fprintf(outfile,"H_update:     None\n");
    else if (optinfo.H_update == OPTInfo::BFGS)
      fprintf(outfile,"H_update:     BFGS\n");
    else if (optinfo.H_update == OPTInfo::MS)
      fprintf(outfile,"H_update:     MS\n");
    else if (optinfo.H_update == OPTInfo::POWELL)
      fprintf(outfile,"H_update:     POWELL\n");
    else if (optinfo.H_update == OPTInfo::BOFILL)
      fprintf(outfile,"H_update:     Bofill\n");

    if (optinfo.empirical_H == OPTInfo::SCHLEGEL)
      fprintf(outfile,"empirical_H:  Schlegel\n");
    else if (optinfo.empirical_H == OPTInfo::FISCHER)
      fprintf(outfile,"empirical_H:  Fischer\n");

    fprintf(outfile,"Max. consecutive line searches: %d\n", optinfo.max_consecutive_line_searches);

    fprintf(outfile,"H_update_use_last: %d\n",optinfo.H_update_use_last);
    fprintf(outfile,"step_limit_cart: %f\n",optinfo.step_limit_cart);
    fprintf(outfile,"step_energy_limit : %f\n",optinfo.step_energy_limit);
    fprintf(outfile,"step_energy_limit_back: %f\n",optinfo.step_energy_limit_back);
    fprintf(outfile,"step_limit: %f\n",optinfo.step_limit);
    fprintf(outfile,"mix_types:     %d\n",optinfo.mix_types);
    fprintf(outfile,"delocalize:    %d\n",optinfo.delocalize);

    if (optinfo.opt_conv == OPTInfo::LOOSE)
      fprintf(outfile,"opt_conv:        loose\n");
    else if (optinfo.opt_conv == OPTInfo::NORMAL)
      fprintf(outfile,"opt_conv:        normal\n");
    else if (optinfo.opt_conv == OPTInfo::TIGHT)
      fprintf(outfile,"opt_conv:        tight\n");
    else if (optinfo.opt_conv == OPTInfo::VERY_TIGHT)
      fprintf(outfile,"opt_conv:        very_tight\n");
    else if (optinfo.opt_conv == OPTInfo::BAKER)
      fprintf(outfile,"opt_conv:        baker\n");
    else if (optinfo.opt_conv == OPTInfo::QCHEM)
      fprintf(outfile,"opt_conv:        qchem\n");
    else if (optinfo.opt_conv == OPTInfo::G03_NORMAL)
      fprintf(outfile,"opt_conv:        g03_normal\n");
    else if (optinfo.opt_conv == OPTInfo::G03_TIGHT)
      fprintf(outfile,"opt_conv:        g03_tight\n");
    else if (optinfo.opt_conv == OPTInfo::G03_VERY_TIGHT)
      fprintf(outfile,"opt_conv:        g03_very_tight\n");
    else if (optinfo.opt_conv == OPTInfo::GENERAL)
      fprintf(outfile,"opt_conv:        general\n");

    fprintf(outfile,"conv_max_force:  %.1e\n",optinfo.conv_max_force);
    fprintf(outfile,"conv_rms_force:  %.1e\n",optinfo.conv_rms_force);
    fprintf(outfile,"conv_max_DE:     %.1e\n",optinfo.conv_max_DE);
    fprintf(outfile,"conv_max_disp:   %.1e\n",optinfo.conv_max_disp);
    fprintf(outfile,"conv_rms_disp:   %.1e\n",optinfo.conv_rms_disp);

    fprintf(outfile,"dertype:       %d\n",optinfo.dertype);
    fprintf(outfile,"numerical dertype: %d\n",optinfo.numerical_dertype);
    fprintf(outfile,"iteration:       %d\n",optinfo.iteration);
    fprintf(outfile,"micro_iteration: %d\n",optinfo.micro_iteration);
    fprintf(outfile,"scale_connectivity: %.3lf\n",optinfo.scale_connectivity);
    fprintf(outfile,"disp_size: %.4lf\n",optinfo.disp_size);
    fprintf(outfile,"bt_max_iter:   %d\n",optinfo.bt_max_iter);
    fprintf(outfile,"bt_max_dq_conv:    %.1e\n",optinfo.bt_dq_conv);
    fprintf(outfile,"bt_max_dx_conv:    %.1e\n",optinfo.bt_dx_conv);
    fprintf(outfile,"fix_interfragment:   %d\n",optinfo.fix_interfragment);
    fprintf(outfile,"fix_intrafragment::  %d\n",optinfo.fix_intrafragment);
    fprintf(outfile,"analytic_interfragment:   %d\n",optinfo.analytic_interfragment);
    fprintf(outfile,"cos_tors_near_1_tol:     %10.6e\n",
        optinfo.cos_tors_near_1_tol);
    fprintf(outfile,"cos_tors_near_neg1_tol: %10.6e\n",
        optinfo.cos_tors_near_neg1_tol);
    fprintf(outfile,"sin_phi_denominator_tol: %10.6e\n",
        optinfo.sin_phi_denominator_tol);
    fflush(outfile);
  }
}

/* POWER raises number to a power */
double power(double x, int y) {
  double tval = 1.0;
  int invert = 0;

  if (y < 0) {
    invert = 1;
    y = abs(y);
  }
  for ( ; y>0; --y) 
    tval *= x;
  if (invert) tval = 1.0/tval;
  return tval;
}

}//} /* namespace psi::optking */

