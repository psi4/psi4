/*! \file
    \ingroup MP2
    \brief Enter brief description of file here 
*/

#include <cstdio>
#include <cstdlib>
#include <libciomr/libciomr.h>
#include <libchkpt/chkpt.h>
#include <libpsio/psio.h>
#include <libqt/qt.h>
#define EXTERN
#include "globals.h"

namespace psi{ namespace mp2{

void get_moinfo(void)
{
  int i;
  
  chkpt_init(PSIO_OPEN_OLD);

  mo.nmo = chkpt_rd_nmo();
  mo.nso = chkpt_rd_nso();
  mo.nao = chkpt_rd_nao();
  mo.nirreps = chkpt_rd_nirreps();
  mo.irreplabels = chkpt_rd_irr_labs();
  mo.mopi = chkpt_rd_orbspi();
  mo.doccpi = chkpt_rd_clsdpi();
  mo.soccpi = chkpt_rd_openpi();
  mo.Enuc = chkpt_rd_enuc();
  mo.Escf = chkpt_rd_escf();

  chkpt_close();
  
  psio_read_entry(CC_INFO,"Reference Wavefunction",(char*)&(params.ref),sizeof(int));
  
  mo.fzdoccpi = init_int_array(mo.nirreps);
  mo.fzvirtpi = init_int_array(mo.nirreps);

  psio_read_entry(CC_INFO,"Frozen Core Orbs Per Irrep",(char*)mo.fzdoccpi,sizeof(int)*mo.nirreps);
  psio_read_entry(CC_INFO,"Frozen Virt Orbs Per Irrep",(char*)mo.fzvirtpi,sizeof(int)*mo.nirreps);
  psio_read_entry(CC_INFO,"No. of Active Orbitals",(char*)&(mo.nactmo),sizeof(int));
  
  if(params.ref == 2) { /** UHF **/

    mo.aoccpi = init_int_array(mo.nirreps);
    mo.boccpi = init_int_array(mo.nirreps);
    mo.avirpi = init_int_array(mo.nirreps);
    mo.bvirpi = init_int_array(mo.nirreps);

    psio_read_entry(CC_INFO,"Active Alpha Occ Orbs Per Irrep",(char*)mo.aoccpi,sizeof(int)*mo.nirreps);
    psio_read_entry(CC_INFO,"Active Beta Occ Orbs Per Irrep",(char*)mo.boccpi,sizeof(int)*mo.nirreps);
    psio_read_entry(CC_INFO,"Active Alpha Virt Orbs Per Irrep",(char*)mo.avirpi,sizeof(int)*mo.nirreps);
    psio_read_entry(CC_INFO,"Active Beta Virt Orbs Per Irrep",(char*)mo.bvirpi,sizeof(int)*mo.nirreps);

    mo.aocc_sym = init_int_array(mo.nactmo);
    mo.bocc_sym = init_int_array(mo.nactmo);
    mo.avir_sym = init_int_array(mo.nactmo);
    mo.bvir_sym = init_int_array(mo.nactmo);

    psio_read_entry(CC_INFO,"Active Alpha Occ Orb Symmetry",(char*)mo.aocc_sym,sizeof(int)*mo.nactmo);
    psio_read_entry(CC_INFO,"Active Beta Occ Orb Symmetry",(char*)mo.bocc_sym,sizeof(int)*mo.nactmo);
    psio_read_entry(CC_INFO,"Active Alpha Virt Orb Symmetry",(char*)mo.avir_sym,sizeof(int)*mo.nactmo);
    psio_read_entry(CC_INFO,"Active Beta Virt Orb Symmetry",(char*)mo.bvir_sym,sizeof(int)*mo.nactmo);

    mo.aocc_off = init_int_array(mo.nirreps);
    mo.bocc_off = init_int_array(mo.nirreps);
    mo.avir_off = init_int_array(mo.nirreps);
    mo.bvir_off = init_int_array(mo.nirreps);

    psio_read_entry(CC_INFO,"Active Alpha Occ Orb Offsets",(char*)mo.aocc_off,sizeof(int)*mo.nirreps);
    psio_read_entry(CC_INFO,"Active Beta Occ Orb Offsets",(char*)mo.bocc_off,sizeof(int)*mo.nirreps);
    psio_read_entry(CC_INFO,"Active Alpha Virt Orb Offsets",(char*)mo.avir_off,sizeof(int)*mo.nirreps);
    psio_read_entry(CC_INFO,"Active Beta Virt Orb Offsets",(char*)mo.bvir_off,sizeof(int)*mo.nirreps);
  
    mo.qt_aocc = init_int_array(mo.nactmo);
    mo.qt_bocc = init_int_array(mo.nactmo);
    mo.qt_avir = init_int_array(mo.nactmo);
    mo.qt_bvir = init_int_array(mo.nactmo);

    psio_read_entry(CC_INFO,"CC->QT Alpha Active Occ Order",(char*)mo.qt_aocc,sizeof(int)*mo.nactmo);
    psio_read_entry(CC_INFO,"CC->QT Beta Active Occ Order",(char*)mo.qt_bocc,sizeof(int)*mo.nactmo);
    psio_read_entry(CC_INFO,"CC->QT Alpha Active Virt Order",(char*)mo.qt_avir,sizeof(int)*mo.nactmo);
    psio_read_entry(CC_INFO,"CC->QT Beta Active Virt Order",(char*)mo.qt_bvir,sizeof(int)*mo.nactmo);
  }
  else { /** RHF or ROHF **/

    mo.occpi = init_int_array(mo.nirreps);
    mo.virpi = init_int_array(mo.nirreps);
    psio_read_entry(CC_INFO,"Active Occ Orbs Per Irrep",(char*)mo.occpi,sizeof(int)*mo.nirreps);
    psio_read_entry(CC_INFO,"Active Virt Orbs Per Irrep",(char*)mo.virpi,sizeof(int)*mo.nirreps);
  
    mo.occ_sym = init_int_array(mo.nactmo);
    mo.vir_sym = init_int_array(mo.nactmo);
    psio_read_entry(CC_INFO,"Active Occ Orb Symmetry",(char*)mo.occ_sym,sizeof(int)*mo.nactmo);
    psio_read_entry(CC_INFO,"Active Virt Orb Symmetry",(char*)mo.vir_sym,sizeof(int)*mo.nactmo);

    mo.occ_off = init_int_array(mo.nirreps);
    mo.vir_off = init_int_array(mo.nirreps);
    psio_read_entry(CC_INFO,"Active Occ Orb Offsets",(char*)mo.occ_off,sizeof(int)*mo.nirreps);
    psio_read_entry(CC_INFO,"Active Virt Orb Offsets",(char*)mo.vir_off,sizeof(int)*mo.nirreps);

    mo.qt_occ = init_int_array(mo.nactmo);
    mo.qt_vir = init_int_array(mo.nactmo);

    psio_read_entry(CC_INFO,"CC->QT Active Occ Order",(char*)mo.qt_occ,sizeof(int)*mo.nactmo);
    psio_read_entry(CC_INFO,"CC->QT Active Virt Order",(char*)mo.qt_vir,sizeof(int)*mo.nactmo);
  }
	      
  mo.virtpi = init_int_array(mo.nirreps);
  for(i=0; i<mo.nirreps; i++)
    mo.virtpi[i] = mo.mopi[i] - mo.doccpi[i] - mo.soccpi[i];

  mo.ndocc = mo.nsocc = mo.nvirt = 0;
  for(i=0; i<mo.nirreps; i++) {
    mo.nfzdocc += mo.fzdoccpi[i];
    mo.nfzvirt += mo.fzvirtpi[i];
    mo.ndocc += mo.doccpi[i];
    mo.nsocc += mo.soccpi[i];
    mo.nvirt += mo.virtpi[i];
  }
  
  /*fprintf(outfile,"\n");
  fprintf(outfile,"\tChkpt Parameters:\n");
  fprintf(outfile,"\t--------------------\n");
  fprintf(outfile,"\tNumber of irreps     = %d\n",mo.nirreps);
  fprintf(outfile,"\tNumber of MOs        = %d\n",mo.nmo);
  fprintf(outfile,"\n");
  fprintf(outfile,
    "\tLabel\tFZDC\tACTD\tDOCC\tSOCC\tACTV\tFZVI\tVIRT\tMOs\n");
  fprintf(outfile,
    "\t-----\t----\t----\t----\t----\t----\t----\t----\t---\n");
  for(i=0; i<mo.nirreps; i++) {
    fprintf(outfile,
    "\t  %s \t  %d\t  %d\t  %d\t  %d\t  %d\t  %d\t  %d\t %d\n",
	    mo.irreplabels[i],mo.fzdoccpi[i],mo.occpi[i],mo.doccpi[i],mo.soccpi[i],
	    mo.virpi[i],mo.fzvirtpi[i],mo.virtpi[i],mo.mopi[i]);
  }*/
  
  fprintf(outfile,"\n");
  fprintf(outfile,"\tNuclear rep. energy     = %20.15f\n",mo.Enuc);
  fprintf(outfile,"\tSCF energy              = %20.15f\n",mo.Escf);
  fflush(outfile);
}

}} /* End namespace */
