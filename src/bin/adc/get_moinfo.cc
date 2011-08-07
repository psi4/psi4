
/*
 *  MOInfo.cc
 *  
 *
 *  Created by M.Saitow on 11/07/14.
 *  Copyright 2010 M.Saitow. All rights reserved.
 *
 */

#include <cstdio>
#include <cstdlib>
#include <libipv1/ip_lib.h>
#include <libciomr/libciomr.h>
#include <libchkpt/chkpt.h>
#include <libchkpt/chkpt.hpp>
#include <libpsio/psio.h>
#include <libqt/qt.h>
#include "MOInfo.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace adc {
	
void get_moinfo(void)
{
	int i, j;
	int isym, sum = 0, irrep;
	int nactive, nirreps;
	chkpt_init(PSIO_OPEN_OLD);
	moinfo.nirreps = chkpt_rd_nirreps();
	moinfo.nmo    = chkpt_rd_nmo();
	moinfo.nso    = chkpt_rd_nso();
	moinfo.nao    = chkpt_rd_nao();
	moinfo.labels = chkpt_rd_irr_labs();
	moinfo.enuc   = chkpt_rd_enuc();
	moinfo.escf   = chkpt_rd_escf();
	moinfo.orbspi = chkpt_rd_orbspi();
	moinfo.clsdpi = chkpt_rd_clsdpi();
	moinfo.openpi = chkpt_rd_openpi();
	//moinfo.usotao = chkpt_rd_usotao();
	moinfo.scf    = chkpt_rd_scf();
	//moinfo.pg     = chkpt_rd_sym_label();
	chkpt_close();
		
	nirreps = moinfo.nirreps;
		
	moinfo.frdocc = init_int_array(nirreps);
	moinfo.frvir   = init_int_array(nirreps);
	psio_read_entry(CC_INFO, "Frozen Core Orbs Per Irrep", (char *) moinfo.frdocc, sizeof(int)*nirreps);
	psio_read_entry(CC_INFO, "Frozen Virt Orbs Per Irrep", (char *) moinfo.frvir, sizeof(int)*nirreps); 
	psio_read_entry(CC_INFO, "No. of Active Orbitals", (char *) &(nactive), sizeof(int));
		
	moinfo.occpi = init_int_array(nirreps);
	moinfo.virpi = init_int_array(nirreps);
	psio_read_entry(CC_INFO, "Active Occ Orbs Per Irrep", (char *)moinfo.occpi, sizeof(int)*nirreps);
	psio_read_entry(CC_INFO, "Active Virt Orbs Per Irrep", (char *)moinfo.virpi, sizeof(int)*nirreps);
      
	moinfo.occ_sym = init_int_array(nactive);
	moinfo.vir_sym = init_int_array(nactive);
	psio_read_entry(CC_INFO, "Active Occ Orb Symmetry", (char *)moinfo.occ_sym, sizeof(int)*nactive);
	psio_read_entry(CC_INFO, "Active Virt Orb Symmetry", (char *)moinfo.vir_sym, sizeof(int)*nactive);
      
	moinfo.occ_off = init_int_array(moinfo.nirreps);
	moinfo.vir_off = init_int_array(moinfo.nirreps);
	psio_read_entry(CC_INFO, "Active Occ Orb Offsets", (char *) moinfo.occ_off, sizeof(int)*moinfo.nirreps);
	psio_read_entry(CC_INFO, "Active Virt Orb Offsets", (char *) moinfo.vir_off, sizeof(int)*moinfo.nirreps);
      
	psio_read_entry(CC_INFO, "Reference Energy", (char *) &(moinfo.eref), sizeof(double));
      
	moinfo.qt_occ = init_int_array(nactive);
	moinfo.qt_vir = init_int_array(nactive);
	psio_read_entry(CC_INFO, "CC->QT Active Occ Order", (char *) moinfo.qt_occ, sizeof(int)*nactive);
	psio_read_entry(CC_INFO, "CC->QT Active Virt Order", (char *) moinfo.qt_vir, sizeof(int)*nactive);	
      
	moinfo.pitzer2qt = init_int_array(moinfo.nmo);
	moinfo.qt2pitzer = init_int_array(moinfo.nmo);
	reorder_qt(moinfo.clsdpi, moinfo.openpi, moinfo.frdocc, moinfo.frvir, 
	moinfo.pitzer2qt, moinfo.orbspi, moinfo.nirreps);
	for(i = 0; i < moinfo.nmo; i++) {
		j = moinfo.pitzer2qt[i];
		moinfo.qt2pitzer[j] = i;
	}	
      
	moinfo.nocc = 0;
	moinfo.nvir = 0;
	for(irrep = 0;irrep < nirreps;irrep++){
		moinfo.nocc += moinfo.occpi[irrep];
		moinfo.nvir += moinfo.virpi[irrep];
	}
      
	moinfo.totdim = 0;
	moinfo.nexc = init_int_array(nirreps);
	for(irrep = 0;irrep < nirreps;irrep++){
		for(isym = 0;isym < nirreps;isym++){
			sum += moinfo.occpi[isym] * moinfo.virpi[isym^irrep];
		}
		moinfo.nexc[irrep] = sum;
		moinfo.totdim += sum;
		sum = 0;
	}
}
    
void cleanup(void)
{
	int i;
	free(moinfo.nexc);
	free(moinfo.orbspi);
	free(moinfo.orbsym);
	free(moinfo.clsdpi);
	free(moinfo.occpi);
	free(moinfo.occ_sym);
	free(moinfo.virpi);
	free(moinfo.vir_sym);
	free(moinfo.evals);
	free(moinfo.qt_occ);
	free(moinfo.frdocc);
	free(moinfo.frvir);
	free(moinfo.openpi);
	free(moinfo.occ_off);
	free(moinfo.vir_off);
	free(moinfo.qt_vir);
	free(moinfo.pitzer2qt);
	free(moinfo.qt2pitzer);
	free(moinfo.pg);
	free_block(moinfo.scf);
	free_block(moinfo.usotao);
	for(i = 0;i < moinfo.nirreps;i++)
		free(moinfo.labels[i]);
}
    
}}


