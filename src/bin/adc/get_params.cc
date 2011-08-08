
/*
 *  get_params.cc
 *  
 *
 *  Created by M.Saitow on 11/07/14.
 *  Copyright 2010 M.Saitow. All rights reserved.
 *
 */

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
//#include <libipv1/ip_lib.h>
#include <libciomr/libciomr.h>
#include <psifiles.h>
#include <libpsio/psio.h>
#include <libchkpt/chkpt.h>
//#include "MOInfo.h"
//#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi {namespace adc {

void get_params(Options &options)
{
	int *mu_irreps;
	int i, errcod, iconv, icrit;
	std::string junk;
    /*
	params.accelerate  = 1;
	params.convergence = 1e-7;
	params.criterion   = 1e-7;
	params.maxiter     = 15;
	params.dim         = 1;
	params.width       = 9;
	params.num_amps    = 5;
	params.order       = 2;
	params.spin        = 0;
    */  
	//fndcor(&(params.memory), infile, outfile);
	params.memory = Process::environment.get_memory();
	guess = (double **)malloc(moinfo.nirreps*sizeof(double*));
	params.d_poles = (double **) malloc(moinfo.nirreps*sizeof(double*));
	params.ocss    = (double **) malloc(moinfo.nirreps*sizeof(double*));
	params.d_ocss  = (double **) malloc(moinfo.nirreps*sizeof(double*));
	params.mu_sqs  = (double ***)malloc(moinfo.nirreps*sizeof(double**));
	params.d_mu_sqs= (double ***)malloc(moinfo.nirreps*sizeof(double**));
	params.rfs     = (double **) malloc(moinfo.nirreps*sizeof(double*));
	params.iter    = (int **)malloc(moinfo.nirreps*sizeof(int*));

	params.accelerate = options.get_int("ACCELERATE");
	params.maxiter    = options.get_int("MAXITER");
	params.dim        = options.get_int("DIM");
	params.width      = options.get_int("WIDTH");
	//iconv             = options.get_int("CONVERGENCE");
	//printf("%d\n", iconv);
	//params.convergence=1.0 * pow(10, -(double)iconv);
	//printf("%lf\n", params.convergence);
	//params.convergence= options.get_double("CONVERGENCE");
	params.num_amps   = options.get_int("NUM_AMPS");

	if(params.num_amps < 0)
		throw PsiException("Num_amps must be greater than 0.\n", __FILE__, __LINE__);
      
	//errcod = ip_data("CRITERION", "%d", &(params.criterion), 0);
	icrit = options.get_int("CRITERION");
	params.criterion = 1.0 * pow(10, -(double)icrit);
	//params.criterion = options.get_double("CRITERION");
      
	params.order = options.get_int("ORDER");
	if(params.order!=1 && params.order!=2)
		throw PsiException("Order of response matrix have to be 1 or 2.\n", __FILE__, __LINE__);
      
	if(options["STATES_PER_IRREP"].size() > 0) {
		i = options["STATES_PER_IRREP"].size();
		if(i != moinfo.nirreps) {
			fprintf(outfile, "Dim. of states_per_irrep vector must be %d\n", moinfo.nirreps);
			throw PsiException("adc input comparison error STATES_PER_IRREP and moinfo.nirreps", __FILE__, __LINE__);
		}
		params.ppi = options.get_int_array("STATES_PER_IRREP");
	}
	else {
		params.ppi = init_int_array(moinfo.nirreps);
		for(i = 0;i < moinfo.nirreps;i++)
			params.ppi[i] = 1;
	}
		
	params.init_dim = init_int_array(moinfo.nirreps);
	for(i = 0;i < moinfo.nirreps;i++){
		if(moinfo.nexc[i] > params.dim * params.ppi[i]) 
			params.init_dim[i] = params.dim * params.ppi[i];
		else 
			params.init_dim[i] = params.ppi[i];
	}
      
	params.start = init_int_array(moinfo.nirreps);
	if(options["START"].size() > 0) {
		i = options["START"].size();
		if (i != moinfo.nirreps) {
			fprintf(outfile, "Dim. of start vector must be %d\n", moinfo.nirreps);
			throw PsiException("adc input comparison error START and moinfo.nirreps", __FILE__, __LINE__);
		}
		params.start = options.get_int_array("START");
	}
	else {
		for(i = 0;i < moinfo.nirreps;i++)
			params.start[i] = 0;
	}

	params.final = init_int_array(moinfo.nirreps);
	if(options["FINAL"].size() > 0){
		i = options["FINAL"].size();
		if(i != moinfo.nirreps){
			fprintf(outfile, "Dim. of FINAL vecor must be %d\n", moinfo.nirreps);
			throw PsiException("adc input comparison error FINAL and moinfo.nirreps", __FILE__, __LINE__);
		}
		params.final = options.get_int_array("FINAL");
	}
	else {
		params.final = init_int_array(moinfo.nirreps);
		for(i = 0;i < moinfo.nirreps;i++)
			params.final[i] = params.ppi[i];
	}
	
	params.spin = options.get_int("SPIN");
	if (params.spin == 0) params.mode = "SINGLET";
	else if (params.spin == 1) params.mode = "TRIPLET";
	else {
		fprintf(stderr, "Spin has to be specified so as to singlet or triplet.\n");
		throw PsiException("Spin state specified is not identifiable.\n", __FILE__, __LINE__);
	}
	
	params.diag_method = "DAVIDSON";
	params.diag_method = options.get_str("DIAG_METHOD");
	if (params.diag_method != "DAVIDSON" && params.diag_method != "EISPACK") {
		fprintf(outfile, "Diag_method should be specified correctly.\n");
		throw PsiException("Only DAVIDSON and EISPACK are available for now.\n", __FILE__, __LINE__);
	}
	
	params.corr_type = "MP2";
	params.corr_type = options.get_str("CORR_TYPE");
	if (params.corr_type != "MP2" && params.corr_type != "CC" && params.corr_type != "PR"){
		fprintf(outfile, "Corr_type should be specified correctly.\n");
		throw PsiException("Only MP2, PR and CC are availeble for the correlated wave function.\n", __FILE__, __LINE__);
	}
	
	junk = options.get_str("REFERENCE");
	if (junk == "RHF") params.ref = 0;
	else if (junk == "ROHF") params.ref = 1;
	else if (junk == "UHF") params.ref = 2;

	params.wfn = options.get_str("WFN");

	params.zero_or_infty = options.get_str("QDPT");
	if (params.zero_or_infty != "INFTY" && params.zero_or_infty != "ZERO") {
		fprintf(stderr, "QDPT has to be specified either ZERO or INFTY.\n");
		throw PsiException("Inconsisitent treatment of quasi-degeneracy.\n", __FILE__, __LINE__);
	}
	      
	mu_irreps = init_int_array(3);
	//mu_irreps = options.get_int_array("MU_IRREPS");
	//moinfo.irrep_x = mu_irreps[0];
	//moinfo.irrep_y = mu_irreps[1];
	//moinfo.irrep_z = mu_irreps[2];
/*	
	if(!strcmp(moinfo.pg, "D2h")){
		moinfo.irrep_x = 7;
		moinfo.irrep_y = 6;
		moinfo.irrep_z = 5;
	}
	else if(!strcmp(moinfo.pg, "D2")){
		moinfo.irrep_x = 3;
		moinfo.irrep_y = 2;
		moinfo.irrep_z = 1;
	}
	else if(!strcmp(moinfo.pg, "C2v")){
		moinfo.irrep_x = 2;
		moinfo.irrep_y = 3;
		moinfo.irrep_z = 0;
	}
	else if(!strcmp(moinfo.pg, "C2h")){
		moinfo.irrep_x = 3;
		moinfo.irrep_y = 3;
		moinfo.irrep_z = 2;
	}
	else if(!strcmp(moinfo.pg, "Cs")){
		moinfo.irrep_x = 0;
		moinfo.irrep_y = 0;
		moinfo.irrep_z = 1;
	}
	else if(!strcmp(moinfo.pg, "C2")){
		moinfo.irrep_x = 1;
		moinfo.irrep_y = 1;
		moinfo.irrep_z = 0;
	}
	else if(!strcmp(moinfo.pg, "Ci")){
		moinfo.irrep_x = 1;
		moinfo.irrep_y = 1;
		moinfo.irrep_z = 1;
	}
	else if(!strcmp(moinfo.pg, "C1 ")){
		moinfo.irrep_x = 0;
		moinfo.irrep_y = 0;
		moinfo.irrep_z = 0;
	}
	else{
		fprintf(outfile, "Point group can not be detected!\n");
		throw PsiException("The point-group can not be detecled.\n", __FILE__, __LINE__);
	}
 */
}

}}

