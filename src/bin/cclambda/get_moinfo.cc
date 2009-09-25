/*! \file
    \ingroup CCLAMBDA
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include <libipv1/ip_lib.h>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>
#include <libchkpt/chkpt.h>
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace cclambda {

/*
** get_moinfo():  Routine to obtain basic orbital information from
** CHKPT and CC_INFO.
**
** T. Daniel Crawford, October 1996.
** Modified by TDC, March 1999.
** Modified for UHF references by TDC, June 2002.
*/

void get_moinfo(void)
{
  int i,j, h, p, q, errcod, nactive, nirreps, sym;
  double ***C, ***Ca, ***Cb;
  psio_address next;

  chkpt_init(PSIO_OPEN_OLD);
  moinfo.nirreps = chkpt_rd_nirreps();
  moinfo.nmo = chkpt_rd_nmo();
  moinfo.nso = chkpt_rd_nso();
  moinfo.nao = chkpt_rd_nao();
  moinfo.iopen = chkpt_rd_iopen();
  moinfo.labels = chkpt_rd_irr_labs();
  moinfo.enuc = chkpt_rd_enuc();
  moinfo.escf = chkpt_rd_escf();
  moinfo.sopi = chkpt_rd_sopi();
  moinfo.orbspi = chkpt_rd_orbspi();
  moinfo.clsdpi = chkpt_rd_clsdpi();
  moinfo.openpi = chkpt_rd_openpi();
  moinfo.phase = chkpt_rd_phase_check();
  chkpt_close();

  sym = 0;
  for (i=0;i<moinfo.nirreps;++i)
    for (j=0;j<moinfo.openpi[i];++j)
      sym = sym ^ i;
  moinfo.sym = sym;

  psio_read_entry(CC_INFO, "Reference Wavefunction", (char *) &(params.ref),
                  sizeof(int));

  /* Get frozen and active orbital lookups from CC_INFO */
  moinfo.frdocc = init_int_array(moinfo.nirreps);
  moinfo.fruocc = init_int_array(moinfo.nirreps);
  psio_read_entry(CC_INFO, "Frozen Core Orbs Per Irrep",
		  (char *) moinfo.frdocc, sizeof(int)*moinfo.nirreps);
  psio_read_entry(CC_INFO, "Frozen Virt Orbs Per Irrep",
		  (char *) moinfo.fruocc, sizeof(int)*moinfo.nirreps);
  
  psio_read_entry(CC_INFO, "No. of Active Orbitals", (char *) &(nactive),
		  sizeof(int)); 

  if(params.ref == 0 || params.ref == 1) { /** RHF or ROHF **/

    moinfo.occpi = init_int_array(moinfo.nirreps);
    moinfo.virtpi = init_int_array(moinfo.nirreps);
    psio_read_entry(CC_INFO, "Active Occ Orbs Per Irrep",
		    (char *) moinfo.occpi, sizeof(int)*moinfo.nirreps);
    psio_read_entry(CC_INFO, "Active Virt Orbs Per Irrep",
		    (char *) moinfo.virtpi, sizeof(int)*moinfo.nirreps);

    moinfo.occ_sym = init_int_array(nactive);
    moinfo.vir_sym = init_int_array(nactive);
    psio_read_entry(CC_INFO, "Active Occ Orb Symmetry",
		    (char *) moinfo.occ_sym, sizeof(int)*nactive);
    psio_read_entry(CC_INFO, "Active Virt Orb Symmetry",
		    (char *) moinfo.vir_sym, sizeof(int)*nactive);

    moinfo.occ_off = init_int_array(moinfo.nirreps);
    moinfo.vir_off = init_int_array(moinfo.nirreps);
    psio_read_entry(CC_INFO, "Active Occ Orb Offsets",
		    (char *) moinfo.occ_off, sizeof(int)*moinfo.nirreps);
    psio_read_entry(CC_INFO, "Active Virt Orb Offsets",
		    (char *) moinfo.vir_off, sizeof(int)*moinfo.nirreps);
  }
  else if(params.ref == 2) { /** UHF **/

    moinfo.aoccpi = init_int_array(moinfo.nirreps);
    moinfo.boccpi = init_int_array(moinfo.nirreps);
    moinfo.avirtpi = init_int_array(moinfo.nirreps);
    moinfo.bvirtpi = init_int_array(moinfo.nirreps);

    psio_read_entry(CC_INFO, "Active Alpha Occ Orbs Per Irrep",
		    (char *) moinfo.aoccpi, sizeof(int)*moinfo.nirreps);
    psio_read_entry(CC_INFO, "Active Beta Occ Orbs Per Irrep",
		    (char *) moinfo.boccpi, sizeof(int)*moinfo.nirreps);
    psio_read_entry(CC_INFO, "Active Alpha Virt Orbs Per Irrep",
		    (char *) moinfo.avirtpi, sizeof(int)*moinfo.nirreps);
    psio_read_entry(CC_INFO, "Active Beta Virt Orbs Per Irrep",
		    (char *) moinfo.bvirtpi, sizeof(int)*moinfo.nirreps);

    moinfo.aocc_sym = init_int_array(nactive);
    moinfo.bocc_sym = init_int_array(nactive);
    moinfo.avir_sym = init_int_array(nactive);
    moinfo.bvir_sym = init_int_array(nactive);

    psio_read_entry(CC_INFO, "Active Alpha Occ Orb Symmetry",
		    (char *) moinfo.aocc_sym, sizeof(int)*nactive);
    psio_read_entry(CC_INFO, "Active Beta Occ Orb Symmetry",
		    (char *) moinfo.bocc_sym, sizeof(int)*nactive);
    psio_read_entry(CC_INFO, "Active Alpha Virt Orb Symmetry",
		    (char *) moinfo.avir_sym, sizeof(int)*nactive);
    psio_read_entry(CC_INFO, "Active Beta Virt Orb Symmetry",
		    (char *) moinfo.bvir_sym, sizeof(int)*nactive);

    moinfo.aocc_off = init_int_array(moinfo.nirreps);
    moinfo.bocc_off = init_int_array(moinfo.nirreps);
    moinfo.avir_off = init_int_array(moinfo.nirreps);
    moinfo.bvir_off = init_int_array(moinfo.nirreps);
    psio_read_entry(CC_INFO, "Active Alpha Occ Orb Offsets",
		    (char *) moinfo.aocc_off, sizeof(int)*moinfo.nirreps);
    psio_read_entry(CC_INFO, "Active Beta Occ Orb Offsets",
		    (char *) moinfo.bocc_off, sizeof(int)*moinfo.nirreps);
    psio_read_entry(CC_INFO, "Active Alpha Virt Orb Offsets",
		    (char *) moinfo.avir_off, sizeof(int)*moinfo.nirreps);
    psio_read_entry(CC_INFO, "Active Beta Virt Orb Offsets",
		    (char *) moinfo.bvir_off, sizeof(int)*moinfo.nirreps);
  }

  /* Build orbsym array (for AO-basis BT2) */
  moinfo.sosym = init_int_array(moinfo.nso);
  for(h=0,q=0; h < moinfo.nirreps; h++)
    for(p=0; p < moinfo.sopi[h]; p++)
      moinfo.sosym[q++] = h;

  if(params.ref == 0 || params.ref == 1) { /** RHF/ROHF **/

    C = (double ***) malloc(moinfo.nirreps * sizeof(double **));
    next = PSIO_ZERO;
    for(h=0; h < moinfo.nirreps; h++) {
      if(moinfo.sopi[h] && moinfo.virtpi[h]) {
        C[h] = block_matrix(moinfo.sopi[h],moinfo.virtpi[h]);
        psio_read(CC_INFO, "RHF/ROHF Active Virtual Orbitals", (char *) C[h][0],
		          moinfo.sopi[h]*moinfo.virtpi[h]*sizeof(double), next, &next);
      }
    }
    moinfo.C = C;
  }
  else if(params.ref == 2) { /** UHF **/

    Ca = (double ***) malloc(moinfo.nirreps * sizeof(double **));
    next = PSIO_ZERO;
    for(h=0; h < moinfo.nirreps; h++) {
      if(moinfo.sopi[h] && moinfo.avirtpi[h]) {
        Ca[h] = block_matrix(moinfo.sopi[h],moinfo.avirtpi[h]);
        psio_read(CC_INFO, "UHF Active Alpha Virtual Orbs", (char *) Ca[h][0],
                  moinfo.sopi[h]*moinfo.avirtpi[h]*sizeof(double), next, &next);
      }
    }
    moinfo.Ca = Ca;


    Cb = (double ***) malloc(moinfo.nirreps * sizeof(double **));
    next = PSIO_ZERO;
    for(h=0; h < moinfo.nirreps; h++) {
      if(moinfo.sopi[h] && moinfo.bvirtpi[h]) {
        Cb[h] = block_matrix(moinfo.sopi[h],moinfo.bvirtpi[h]);
        psio_read(CC_INFO, "UHF Active Beta Virtual Orbs", (char *) Cb[h][0],
                  moinfo.sopi[h]*moinfo.bvirtpi[h]*sizeof(double), next, &next);
      }
    }
    moinfo.Cb = Cb;
  }

  /* Adjust clsdpi array for frozen orbitals */
  for(i=0; i < moinfo.nirreps; i++)
    moinfo.clsdpi[i] -= moinfo.frdocc[i];

  moinfo.uoccpi = init_int_array(moinfo.nirreps);
  for(i=0; i < moinfo.nirreps; i++)
    moinfo.uoccpi[i] = moinfo.orbspi[i] - moinfo.clsdpi[i] -
      moinfo.openpi[i] - moinfo.fruocc[i] -
      moinfo.frdocc[i];

  if(params.ref == 0) {
    moinfo.nvirt = 0;
    for(h=0; h < moinfo.nirreps; h++) moinfo.nvirt += moinfo.virtpi[h];
  }

  psio_read_entry(CC_INFO, "Reference Energy", (char *) &(moinfo.eref),
		  sizeof(double));

  fprintf(outfile,"\n\tNuclear Rep. energy (chkpt)   = %20.15f\n",moinfo.enuc);
  fprintf(outfile,  "\tReference           (chkpt)   = %d\n",params.ref);
  fprintf(outfile,  "\tSCF energy          (chkpt)   = %20.15f\n",moinfo.escf);
  fprintf(outfile,  "\tReference energy    (CC_INFO) = %20.15f\n",moinfo.eref);

  fflush(outfile);
}

/* Frees memory allocated in get_moinfo() and dumps some info. */
void cleanup(void)
{
  int i, h;

  psio_write_entry(CC_INFO, "Lambda Pseudoenergy", (char *) &(moinfo.lcc),
		   sizeof(double));

  if(params.ref == 0 || params.ref == 1) {
    for(h=0; h < moinfo.nirreps; h++)
      if(moinfo.sopi[h] && moinfo.virtpi[h]) free_block(moinfo.C[h]);
    free(moinfo.C);
  }
  else if(params.ref == 2) {
    for(h=0; h < moinfo.nirreps; h++)
      if(moinfo.sopi[h] && moinfo.avirtpi[h]) free_block(moinfo.Ca[h]);
    free(moinfo.Ca);
    for(h=0; h < moinfo.nirreps; h++)
      if(moinfo.sopi[h] && moinfo.bvirtpi[h]) free_block(moinfo.Cb[h]);
    free(moinfo.Cb);
  }

  free(moinfo.sopi);
  free(moinfo.sosym);
  free(moinfo.orbspi);
  free(moinfo.clsdpi);
  free(moinfo.openpi);
  free(moinfo.uoccpi);
  free(moinfo.fruocc);
  free(moinfo.frdocc);
  for(i=0; i < moinfo.nirreps; i++)
    free(moinfo.labels[i]);
  free(moinfo.labels);
  if(params.ref == 2) {
    free(moinfo.aocc_sym);
    free(moinfo.bocc_sym);
    free(moinfo.avir_sym);
    free(moinfo.bvir_sym);
    free(moinfo.aocc_off);
    free(moinfo.bocc_off);
    free(moinfo.avir_off);
    free(moinfo.bvir_off);
    free(moinfo.aoccpi);
    free(moinfo.boccpi);
    free(moinfo.avirtpi);
    free(moinfo.bvirtpi);
  }
  else {
    free(moinfo.occ_sym);
    free(moinfo.vir_sym);
    free(moinfo.occ_off);
    free(moinfo.vir_off);
    free(moinfo.occpi);
    free(moinfo.virtpi);
  }
}


}} // namespace psi::cclambda
