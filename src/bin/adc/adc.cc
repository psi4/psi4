
/*
 *  adc.cc
 *  
 *
 *  Created by M.Saitow on 11/07/14.
 *  Copyright 2010 M.Saitow. All rights reserved.
 *
 */

#include <cstdlib>
#include <cstring>
#include <cmath>
#include <libqt/qt.h>
#include <libipv1/ip_lib.h>
#include <libpsio/psio.h>
#include <libciomr/libciomr.h>
#include <libchkpt/chkpt.h>
#include <libdpd/dpd.h>
#include <psifiles.h>
#include <physconst.h>
#include "Params.h"
#include "MOInfo.h"
#include "globals.h"

namespace psi { namespace adc {

void init_psi(void);
void init_ints(void);
void exit_psi(void);
void get_moinfo(void);
void get_params(Options &options);
void cleanup(void);
void show_params(void);
void title(void);
void make_diag(void);
void make_pp(void);
void renormalize(void);
void diag_adc(void);
void mp2(void);
void ccamps(void);
void amp_write_T1(dpdfile2*, int, FILE*);
void transpert(void);
void sort_pert(void);
double diagonalize(int, int, double);
double accelerate(int, int);
double calc_ocss(char*, int, int, int, double); 
    
int **cacheprep_rhf(int, int *);
void cachedone_rhf(int **);

#define CUTOFF 1.0e-7

PsiReturnType
adc(Options &options)
{
	using namespace psi::adc;
	
  char lbl[32], *axis;
  int i;
  int **cachelist, *cachefiles;
  int irrep, root, iter, nconv;;
  double omega, omega_diff, epsilon, value, rrf=1;
  dpdfile2 V;
  FILE *iter_adc;

  init_psi();
  title();
  get_moinfo();
  get_params(options);
  cachefiles = init_int_array(PSIO_MAXUNIT);
  cachelist = cacheprep_rhf(0, cachefiles);
  dpd_init(0, moinfo.nirreps, params.memory, 0, cachefiles,
	   cachelist, NULL, 2, moinfo.occpi, moinfo.occ_sym,
	   moinfo.virpi, moinfo.vir_sym);	
  
  timer_init();
  show_params();
  fflush(outfile);
  make_pp(); 
  diag_adc();
  
  if(params.order == 2){
    
    if((params.corr_type == "MP2") || (params.corr_type == "PR"))
      mp2();
    else
      ccamps();
    
    renormalize();
    init_ints();
    make_diag();
    
    axis = (char * )malloc(3*sizeof(char));
    axis[0] = 'X'; axis[1] = 'Y'; axis[2] = 'Z';
    
    fprintf(outfile, "\n");
    fprintf(outfile, "\tRHF-");
    if(params.corr_type != "MP2" && params.corr_type == "PR")
      fprintf(outfile, "PR-");
    else if(params.corr_type == "CC")
      fprintf(outfile, "CC-");
    if(params.zero_or_infty == "INFTY")
      fprintf(outfile, "ADC(2)");
    else
      fprintf(outfile, "CIS(2X)");
    if(params.spin == 0) fprintf(outfile, " Singlet");
    else fprintf(outfile, " Triplet");
    fprintf(outfile, " Excitation Energies:\n");
    
    fprintf(outfile, "\t------------------------------------\n\n");
    fprintf(outfile, "\t  Root Irrep       Hartree          eV          cm-1       rf\n");
    fprintf(outfile, "\t  ---- ----- ------------------  ---------  ----------  ---------\n");
    
    for(irrep = 0;irrep < moinfo.nirreps;irrep++){
      ffile(&iter_adc, "iter.dat", 1);
      fprintf(iter_adc, "\tirrep  %s\n", moinfo.labels[irrep]);
      fclose(iter_adc);
      params.d_poles[irrep] = init_array(params.ppi[irrep]);
      params.rfs[irrep]     = init_array(params.ppi[irrep]);
      params.iter[irrep]    = init_int_array(params.ppi[irrep]);
      for(root = params.start[irrep];root < params.final[irrep];root++){
	
	if(params.zero_or_infty == "ZERO"){
	  omega = diagonalize(irrep, root, 0);
	} else{
	  omega = guess[irrep][root];
	  nconv = 0;
	  for(iter = 1;iter < params.maxiter;iter++){
	    epsilon = diagonalize(irrep, root, omega);
	    ffile(&iter_adc, "iter.dat", 1);
	    if(fabs(epsilon) < CUTOFF) nconv++;
	    if(nconv >= 2) {
	      fprintf(outfile, "\t%dth root in %s is oscillating...\n", root, moinfo.labels[irrep]);
	      fprintf(iter_adc, "\t%dth root in %s is oscillating...\n", root, moinfo.labels[irrep]);
	      iter = -1;
	      break;
	    }					
	    if(params.accelerate)
	      rrf = accelerate(root, irrep);
	    omega_diff = (omega-epsilon) / rrf;
	    fprintf(iter_adc, "\t%d. %10.7lf\n", iter, omega*_hartree2ev);
	    fclose(iter_adc);
	    if(fabs(omega_diff) <= params.criterion) break;
	    else
	      omega -= omega_diff;
	  }
	  if(iter>0){
	    ffile(&iter_adc, "iter.dat", 1);
	    fprintf(iter_adc, "\t%dth pole in %d iteration: %10.7lf %10.7lf r.f. %10.7lf\n", 
		    root+1, iter, omega, omega*_hartree2ev, 1/rrf);
	    fclose(iter_adc);
	  }
	  
	}
	sprintf(lbl, "V^(%d)[%d]", root, irrep);
	params.d_poles[irrep][root] = omega;
	params.rfs[irrep][root] = 1 / rrf;
	params.iter[irrep][root] = iter;
	//params.d_ocss[irrep][root] = calc_ocss(lbl, irrep, 1, root, omega);
	
	value = params.d_poles[irrep][root];
	fprintf(outfile,   "\tState %4d   %3s %18.14f  %9.5f  %10.2f   %10.8f\n", irrep, moinfo.labels[irrep],
		value, value*_hartree2ev, value*_hartree2wavenumbers, params.rfs[irrep][root]);
	fprintf(outfile, "\n\tLargest components of singlet excited wave function #%d/#%d:\n",
		irrep, root);
	dpd_file2_init(&V, CC_OEI, irrep, 0, 1, lbl);
	amp_write_T1(&V, params.num_amps, outfile);
	dpd_file2_close(&V);
	if ((params.iter[irrep][root]<params.maxiter) && (params.iter[irrep][root]>0)){
	  fprintf(outfile, "\tConverged in %2d iter.\n", params.iter[irrep][root]);
	  //fprintf(outfile, "\tComponents of squared dipole moment:\n");
	  //fprintf(outfile, "\t");
	  //for(i = 0;i < 3;i++)
	    //fprintf(outfile, "%c: %10.7f ", axis[i], params.d_mu_sqs[irrep][root][i]);
	  //fprintf(outfile, "\n\tOscillator strength is %10.7f.\n\n", params.d_ocss[irrep][root]);
	  fflush(outfile);					
	}
	else if(params.zero_or_infty == "INFTY")
	  fprintf(outfile, "\tI'm so sorry for N.C.\n");
	fprintf(outfile, "\n");
	
	fflush(outfile);
      }
    }
  }
  timer_done();
  dpd_close(0);
  cachedone_rhf(cachelist);
  free(cachefiles);	
  cleanup();
  exit_psi();

  return Success;
}

	/*
extern "C" {
  const char *gprgid(void) 
  {
    const char *prgid = "ADC";
    return(prgid);
  }
}
*/
	 
void init_psi(void)
{
//      char *progid;      
//     progid = (char *)malloc(strlen(gprgid()+2));
//      sprintf(progid, ":%s", gprgid());
//psi_start(&infile, &outfile, &psi_file_prefix, argc-1, argv+1, 0);
//ip_cwk_add(progid)
//free(progid);
	tstart();
//psio_init(); psio_ipv1_config();
/*      
	if(params.corr_type == "CC");
		psio_open(CC_TAMPS, 1);
	psio_open(CC_INFO, 1);
	psio_open(CC_OEI, 1);
	psio_open(CC_CINTS, 1);
	psio_open(CC_DINTS, 1);
	psio_open(CC_EINTS, 1);
	psio_open(CC_FINTS, 1);
	psio_open(CC_DENOM, 1);
	psio_open(CC_MISC, 1); //0
	psio_open(CC_TMP0, 1); //0
	psio_open(CC_TMP1, 1); //0
*/
	for (int i = CC_MIN;i <= CC_MAX;i++)
		psio_open(i,1);
}
	
void exit_psi(void)
{
	/*
	if(params.corr_type == "CC");
		psio_close(CC_TAMPS, 1);
	psio_close(CC_INFO, 1);
	psio_close(CC_OEI, 1);
	psio_close(CC_CINTS, 1);
	psio_close(CC_DINTS, 1);
	psio_close(CC_EINTS, 1);
	psio_close(CC_FINTS, 1);
	psio_close(CC_DENOM, 1);
	psio_close(CC_MISC, 1);
	psio_close(CC_TMP0, 0);
	psio_close(CC_TMP1, 0);
	*/
	//psio_done();
	tstop();
	//psi_stop(infile, outfile, psi_file_prefix);
}
    
void title(void)
{
	fprintf(outfile, "\n");
	fprintf(outfile, "\t***************************************\n");
	fprintf(outfile, "\t* ADC: ADC code for excitation energy *\n");
	fprintf(outfile, "\t*                                     *\n");
	fprintf(outfile, "\t*           Masaaki Saitow            *\n");
	fprintf(outfile, "\t*              Aug 2011               *\n");
	fprintf(outfile, "\t*              Kanazawa               *\n");
	fprintf(outfile, "\t***************************************\n");
	fprintf(outfile, "\n");
	fflush(outfile);
	
	return;
}
    
void show_params(void)
{
	int i;
	fprintf(outfile, "\tInput parameters:\n");
	fprintf(outfile, "\t--------------------\n");
	fprintf(outfile, "\tWavefunction               = %s\n", params.wfn.c_str());
	fprintf(outfile, "\tLevels of QDPT             = %s\n", params.zero_or_infty.c_str());
	fprintf(outfile, "\tReference WFN              = %s\n",
	(params.ref==0) ? "RHF" : "UHF");
	fprintf(outfile, "\tTotal dimension            = %3d\n", moinfo.totdim);
	fprintf(outfile, "\tInit_dim of Ritz space     = %3d *\n",  params.dim);
	fprintf(outfile, "\t Max_dim of Ritz space     = %3d *\n",  params.dim+params.width);
	fprintf(outfile, "\tOrder of the fluctuation   = %3d\n", params.order);
	fprintf(outfile, "\tNum_amps                   = %3d\n", params.num_amps);
	if(params.order == 2){
		fprintf(outfile, "\tMaximum number of iteration    = %3d\n", params.maxiter);
		fprintf(outfile, "\tConvergence criterion(a.u.)    = %e\n", params.criterion);
	}
	fprintf(outfile, "\tNewton-Raphson acceleration    = %s\n", (params.accelerate) ? "TRUE" : "NONE");
      
	fprintf(outfile, "\tIrrep for each of X, Y, and Z  = [%4s %4s %4s]\n", moinfo.labels[moinfo.irrep_x]
	      , moinfo.labels[moinfo.irrep_y]
	      , moinfo.labels[moinfo.irrep_z]);	
      
	fprintf(outfile, "\tMemory (MiB)                 = %5.1f\n",params.memory/1e6);
	fprintf(outfile, "\tRoots sought per irrep       = ");
	for(i = 0;i < moinfo.nirreps;++i) fprintf(outfile, "%8d", params.ppi[i]);
	fprintf(outfile, "\n");
	fprintf(outfile, "\tStart guess for each irrep   = ");
	for(i = 0;i < moinfo.nirreps;++i) fprintf(outfile, "%8d", params.start[i]);
	fprintf(outfile, "\n");
	fprintf(outfile, "\tFinal guess for each irrep   = ");
	for(i = 0;i < moinfo.nirreps;++i) fprintf(outfile, "%8d", params.final[i]);
	fprintf(outfile, "\n");	
	fprintf(outfile, "\tSub-dimension for each irrep = ");
	for(i = 0;i < moinfo.nirreps;i++) fprintf(outfile, "%8d", moinfo.nexc[i]);
	fprintf(outfile, "\n");
	fprintf(outfile, "\tMultiplicity of excited states    = %s\n", params.mode.c_str());
	fprintf(outfile, "\tCorrelation type                  = %s\n", params.corr_type.c_str());
	fprintf(outfile, "\tDiagonalization method            = %s\n", params.diag_method.c_str());
	fprintf(outfile, "\n");
	fprintf(outfile, "\tNuclear Rep. energy (chkpt)   = %20.15f\n",moinfo.enuc);
	fprintf(outfile, "\tSCF energy          (chkpt)   = %20.15f\n",moinfo.escf);
	fprintf(outfile, "\tReference energy    (file100) = %20.15f\n",moinfo.eref);
	fprintf(outfile, "\n");
}
    
}}

