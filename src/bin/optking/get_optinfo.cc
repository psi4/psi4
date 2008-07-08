/*! \file
    \ingroup OPTKING
    \brief GET_OPTINFO   reads optimization parameters from input.dat
*/

#define EXTERN
#include "optking.h"
#undef EXTERN

namespace psi { namespace optking {

double power(double x, int y);

void get_optinfo() {
  int a, i, cnt, cnt2, nallatom, errcod;
  char *junk;

  optinfo.print_lvl = 1;
  ip_data("PRINT_LVL", "%d", &(optinfo.print_lvl),0);

  optinfo.scale_connectivity = 1.3;
  ip_data("SCALE_CONNECTIVITY", "%f",&(optinfo.scale_connectivity),0);

  chkpt_init(PSIO_OPEN_OLD);
  optinfo.natom = chkpt_rd_natom();
  chkpt_close();

/*
  optinfo.iteration = 0;
  optinfo.micro_iteration = 0;
  open_PSIF();
  if (psio_tocscan(PSIF_OPTKING, "Iteration") != NULL)
    psio_read_entry(PSIF_OPTKING, "Iteration",
        (char *) &(optinfo.iteration),sizeof(int));
  if (psio_tocscan(PSIF_OPTKING, "Micro_iteration") != NULL)
    psio_read_entry(PSIF_OPTKING, "Micro_iteration",
        (char *) &(optinfo.micro_iteration),sizeof(int));
  close_PSIF();

  optinfo.dertype = 0;
  if (ip_exist("DERTYPE",0)) {
    errcod = ip_string("DERTYPE", &(junk),0);
    if (errcod != IPE_OK) optinfo.dertype = 0;
    else if(!strcmp(junk,"NONE")) optinfo.dertype = 0;
    else if(!strcmp(junk,"FIRST")) optinfo.dertype = 1;
    else if(!strcmp(junk,"SECOND")) optinfo.dertype = 2;
    else if(!strcmp(junk,"RESPONSE")) optinfo.dertype = 3;
    else {
      printf("Invalid value of input keyword DERTYPE: %s\n", junk);
      exit(PSI_RETURN_FAILURE);
    }
  }

  errcod = ip_string("WFN", &(optinfo.wfn), 0);
  errcod = ip_string("JOBTYPE", &(optinfo.jobtype), 0);
 
  optinfo.numerical_dertype = 0;
  if (ip_exist("DERTYPE",0)) {
    if ( (strcmp(optinfo.jobtype,"OPT")==0) && (optinfo.dertype==0) ) 
      optinfo.numerical_dertype = 1;
  }

  optinfo.redundant = 1; optinfo.delocalize = 0;
  if ((optinfo.mode == MODE_DISP_IRREP) || (optinfo.mode == MODE_DISP_NOSYMM) ) 
    { optinfo.redundant = 0; optinfo.delocalize =1; }
  if ( optinfo.numerical_dertype == 1 )
    { optinfo.redundant = 0; optinfo.delocalize =1; }
  ip_boolean("DELOCALIZE", &(optinfo.delocalize),0);
  if (optinfo.delocalize)
    optinfo.redundant = 0;

  optinfo.frag_dist_rho = 0;
  ip_boolean("FRAGMENT_DISTANCE_INVERSE", &(optinfo.frag_dist_rho),0);
  optinfo.fix_intrafragment = 0;
  ip_boolean("FIX_INTRAFRAGMENT", &(optinfo.fix_intrafragment),0);
  optinfo.fix_interfragment = 0;
  ip_boolean("FIX_INTERFRAGMENT", &(optinfo.fix_interfragment),0);
  optinfo.freeze_intrafragment = 0;
  ip_boolean("FREEZE_INTRAFRAGMENT", &(optinfo.freeze_intrafragment),0);

  optinfo.print_hessian = 0;
  ip_boolean("PRINT_HESSIAN", &(optinfo.print_hessian),0);
  optinfo.optimize = 1;
  if (ip_exist("DISPLACEMENTS",0)) {
    optinfo.do_only_displacements = 1;
    optinfo.mode = MODE_DISP_USER;
    optinfo.optimize = 0;
  }

  optinfo.freq_irrep = -1;
  if (ip_exist("FREQ_IRREP",0)) {
    ip_data("FREQ_IRREP","%d",&(optinfo.freq_irrep),0);
  }

  optinfo.bfgs = 1;
*/
/*
  ip_boolean("BFGS",&(optinfo.bfgs),0);
  // ACS (11/06) Are we getting our energies from outside of PSI3?
  optinfo.external_energies = 0;
  ip_boolean("EXTERNAL_ENERGIES",&(optinfo.external_energies),0);

  optinfo.bfgs_use_last = 6;
  ip_data("BFGS_USE_LAST","%d",&(optinfo.bfgs_use_last),0);
  optinfo.mix_types = 1;
  ip_boolean("MIX_TYPES",&(optinfo.mix_types),0);

  ip_data("POINTS","%d",&(optinfo.points),0);
*/

  /* takes values of 1,2,3 for x,y,z for location of first dummy of linear bend*/
/*
  optinfo.dummy_axis_1 = 2;
  ip_data("DUMMY_AXIS_1","%d",&(optinfo.dummy_axis_1),0);
  optinfo.dummy_axis_1 -= 1;
  optinfo.dummy_axis_2 = 3;
  ip_data("DUMMY_AXIS_2","%d",&(optinfo.dummy_axis_2),0);
  optinfo.dummy_axis_2 -= 1;

  optinfo.zmat = 0;
  if (ip_exist("ZMAT",0)) optinfo.zmat = 1;

  optinfo.zmat_simples = 0;
  ip_boolean("ZMAT_SIMPLES",&(optinfo.zmat_simples),0);

  optinfo.test_B = 0;
  ip_boolean("TEST_B",&(optinfo.test_B),0);

  a = 5;
  ip_data("CONV","%d",&a,0);
  optinfo.conv = power(10.0, -1*a);

  a= 5;
  ip_data("EV_TOL","%d",&a,0);
  optinfo.ev_tol = power(10.0, -1*a);

  optinfo.scale_connectivity = 1.3;
  ip_data("SCALE_CONNECTIVITY","%lf",&(optinfo.scale_connectivity),0);

  optinfo.disp_size = 0.005;
  ip_data("EDISP","%lf",&(optinfo.disp_size),0);
*/

  /* back-transformation parameters */
/*
  optinfo.bt_max_iter = 60;
  ip_data("BT_MAX_ITER","%d",&(optinfo.bt_max_iter),0);
  a = 2; // some minimal level of this seems necessary for butane
  ip_data("BT_DQ_CONV","%d",&a,0);
  optinfo.bt_dq_conv  = power(10.0, -1*a);
  a = 7;
  ip_data("BT_DX_CONV","%d",&a,0);
  optinfo.bt_dx_conv  = power(10.0, -1*a);
*/

  /* These determing how seriously to take torsions very near 180 */
/*
  optinfo.cos_tors_near_1_tol    =  0.9999999;
  ip_data("COS_TORS_NEAR_1_TOL","%lf",&(optinfo.cos_tors_near_1_tol),0);
  optinfo.cos_tors_near_neg1_tol = -0.9999999;
  ip_data("COS_TORS_NEAR_NEG1_TOL","%lf",&(optinfo.cos_tors_near_neg1_tol),0);
  optinfo.sin_phi_denominator_tol = 1.0e-7;
  ip_data("SIN_PHI_DENOMINATOR_TOL","%lf",&(optinfo.sin_phi_denominator_tol),0);
*/

  /* Compute dummy atom lookup arrays */
/*
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
*/

  /*
  fprintf(outfile,"\nIteration: %d     ",optinfo.iteration+1);
  if (optinfo.numerical_dertype > 0)
    fprintf(outfile,"Micro_iteration: %d",optinfo.micro_iteration+1);
  fprintf(outfile,"\n");
  */
  if (optinfo.print_lvl > 0) {
    fprintf(outfile,"\n+++ Optinfo Parameters +++\n");
    fprintf(outfile,"print_lvl:  %d\n",optinfo.print_lvl);
    fprintf(outfile,"scale_connectivity: %.3lf\n",optinfo.scale_connectivity);
/*
    fprintf(outfile,"optimize:      %d\n",optinfo.optimize);
    fprintf(outfile,"zmat:          %d\n",optinfo.zmat);
    fprintf(outfile,"dummy_axis_1:    %d\n",optinfo.dummy_axis_1);
    fprintf(outfile,"dummy_axis_2:    %d\n",optinfo.dummy_axis_2);
    fprintf(outfile,"sacc fd points:  %d\n",optinfo.points);
    fprintf(outfile,"zmat_simples:  %d\n",optinfo.zmat_simples);
    fprintf(outfile,"redundant:     %d\n",optinfo.redundant);
    fprintf(outfile,"bfgs:          %d\n",optinfo.bfgs);
    fprintf(outfile,"bfgs_use_last: %d\n",optinfo.bfgs_use_last);
    fprintf(outfile,"mix_types:     %d\n",optinfo.mix_types);
    fprintf(outfile,"delocalize:    %d\n",optinfo.delocalize);
    fprintf(outfile,"conv:          %.1e\n",optinfo.conv);
    fprintf(outfile,"dertype:       %d\n",optinfo.dertype);
    fprintf(outfile,"numerical dertype: %d\n",optinfo.numerical_dertype);
    fprintf(outfile,"iteration:       %d\n",optinfo.iteration);
    fprintf(outfile,"micro_iteration: %d\n",optinfo.micro_iteration);
    fprintf(outfile,"disp_size: %.4lf\n",optinfo.disp_size);
    fprintf(outfile,"bt_max_iter:   %d\n",optinfo.bt_max_iter);
    fprintf(outfile,"bt_max_dq_conv:    %.1e\n",optinfo.bt_dq_conv);
    fprintf(outfile,"bt_max_dx_conv:    %.1e\n",optinfo.bt_dx_conv);
    fprintf(outfile,"freeze_intrafragment:    %d\n",optinfo.freeze_intrafragment);
    fprintf(outfile,"cos_tors_near_1_tol:     %10.6e\n",
        optinfo.cos_tors_near_1_tol);
    fprintf(outfile,"cos_tors_near_neg1_tol: %10.6e\n",
        optinfo.cos_tors_near_neg1_tol);
    fprintf(outfile,"sin_phi_denominator_tol: %10.6e\n",
        optinfo.sin_phi_denominator_tol);
*/
    fflush(outfile);
  }
}

// POWER raises number to a power */
/*
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
*/

}} /* namespace psi::optking */

