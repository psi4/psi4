
/*
 *  MOInfo.h
 *  
 *
 *  Created by M.Saitow on 11/07/15.
 *  Copyright 2010 M.Saitow. All rights reserved.
 *
 */

#ifndef _psi_src_bin_adc_MOInfo_h_
#define _psi_src_bin_adc_MOInfo_h_

struct MOInfo {
	int nirreps;
	int nocc;
	int nvir;
	int totdim;
	int nmo;
	int nso;
	int nao;
	int irrep_x;
	int irrep_y;
	int irrep_z;
	int *nexc;
	int *orbspi;        
	int *orbsym;      
	int *clsdpi;         
	int *occpi;
	int *occ_off;
	int *occ_sym;
	int *qt_occ;
	int *frdocc;
	int *openpi;
	int *virpi;
	int *vir_off;
	int *vir_sym;
	int *qt_vir;
	int *frvir;
	int *pitzer2qt;
	int *qt2pitzer;
      
	char *pg;
	char **labels;     
		
	double enuc;
	double escf;
	double eref;
	double *evals;
	double **scf;
	double **usotao;
};

#endif