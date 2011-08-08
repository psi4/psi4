
/*
 *  Params.h
 *  
 *
 *  Created by M.Saitow on 10/07/27.
 *  Copyright 2010 M.Saitow. All rights reserved.
 *
 */

#ifndef _psi_src_bin_adc_Params_h_
#define _psi_src_bin_adc_Params_h_

struct Params {
	int accelerate;
	int order;
	int ref;
	int dim;
	int width;
	int spin;
	int maxiter;
	int num_amps;
	long int memory;
	int *ppi;
	int *init_dim;
	int *start;
	int *final;
	int **iter;
	//double convergence;
	double criterion;
	double **d_poles;
	double **ocss;
	double **d_ocss;
	double **rfs;
	double ***mu_sqs;
	double ***d_mu_sqs;
	std::string wfn;
	std::string diag_method;
	std::string mode;
	std::string zero_or_infty;
	std::string corr_type;
};
   
#endif
