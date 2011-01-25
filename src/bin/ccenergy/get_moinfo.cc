/*! \file
    \ingroup CCENERGY
    \brief Enter brief description of file here
*/
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <libciomr/libciomr.h>
#include <libchkpt/chkpt.h>
#include <libpsio/psio.h>
#include <libqt/qt.h>
#include <psifiles.h>
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccenergy {

/*
** get_moinfo():  Routine to obtain basic orbital information from
** CHKPT and CC_INFO.
**
** T. Daniel Crawford, October 1996
** Modified by TDC, March 1999
*/

void get_moinfo(void)
{
  int i, j, h, p, q, errcod, nactive, nirreps;
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

  nirreps = moinfo.nirreps;

  psio_read_entry(CC_INFO, "Reference Wavefunction", (char *) &(params.ref),
          sizeof(int));

  /* Get frozen and active orbital lookups from CC_INFO */
  moinfo.frdocc = init_int_array(nirreps);
  moinfo.fruocc = init_int_array(nirreps);
  psio_read_entry(CC_INFO, "Frozen Core Orbs Per Irrep",
          (char *) moinfo.frdocc, sizeof(int)*nirreps);
  psio_read_entry(CC_INFO, "Frozen Virt Orbs Per Irrep",
          (char *) moinfo.fruocc, sizeof(int)*nirreps);

  psio_read_entry(CC_INFO, "No. of Active Orbitals", (char *) &(nactive),
          sizeof(int));

  if(params.ref == 2) { /** UHF **/

    moinfo.aoccpi = init_int_array(nirreps);
    moinfo.boccpi = init_int_array(nirreps);
    moinfo.avirtpi = init_int_array(nirreps);
    moinfo.bvirtpi = init_int_array(nirreps);

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

    /* Get CC->QT and QT->CC active occupied and virtual reordering arrays */
    moinfo.qt_aocc = init_int_array(nactive);
    moinfo.qt_bocc = init_int_array(nactive);
    moinfo.qt_avir = init_int_array(nactive);
    moinfo.qt_bvir = init_int_array(nactive);
    moinfo.cc_aocc = init_int_array(nactive);
    moinfo.cc_bocc = init_int_array(nactive);
    moinfo.cc_avir = init_int_array(nactive);
    moinfo.cc_bvir = init_int_array(nactive);

    psio_read_entry(CC_INFO, "CC->QT Alpha Active Occ Order",
            (char *) moinfo.qt_aocc, sizeof(int)*nactive);
    psio_read_entry(CC_INFO, "CC->QT Beta Active Occ Order",
            (char *) moinfo.qt_bocc, sizeof(int)*nactive);
    psio_read_entry(CC_INFO, "CC->QT Alpha Active Virt Order",
            (char *) moinfo.qt_avir, sizeof(int)*nactive);
    psio_read_entry(CC_INFO, "CC->QT Beta Active Virt Order",
            (char *) moinfo.qt_bvir, sizeof(int)*nactive);
    psio_read_entry(CC_INFO, "QT->CC Alpha Active Occ Order",
            (char *) moinfo.cc_aocc, sizeof(int)*nactive);
    psio_read_entry(CC_INFO, "QT->CC Beta Active Occ Order",
            (char *) moinfo.cc_bocc, sizeof(int)*nactive);
    psio_read_entry(CC_INFO, "QT->CC Alpha Active Virt Order",
            (char *) moinfo.cc_avir, sizeof(int)*nactive);
    psio_read_entry(CC_INFO, "QT->CC Beta Active Virt Order",
            (char *) moinfo.cc_bvir, sizeof(int)*nactive);


  }
  else { /** RHF or ROHF **/

    moinfo.occpi = init_int_array(nirreps);
    moinfo.virtpi = init_int_array(nirreps);
    psio_read_entry(CC_INFO, "Active Occ Orbs Per Irrep",
            (char *) moinfo.occpi, sizeof(int)*nirreps);
    psio_read_entry(CC_INFO, "Active Virt Orbs Per Irrep",
            (char *) moinfo.virtpi, sizeof(int)*nirreps);

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

    /* Get CC->QT and QT->CC active occupied and virtual reordering arrays */
    moinfo.qt_occ = init_int_array(nactive);
    moinfo.qt_vir = init_int_array(nactive);
    moinfo.cc_occ = init_int_array(nactive);
    moinfo.cc_vir = init_int_array(nactive);
    psio_read_entry(CC_INFO, "CC->QT Active Occ Order",
            (char *) moinfo.qt_occ, sizeof(int)*nactive);
    psio_read_entry(CC_INFO, "CC->QT Active Virt Order",
            (char *) moinfo.qt_vir, sizeof(int)*nactive);
    psio_read_entry(CC_INFO, "QT->CC Active Occ Order",
            (char *) moinfo.cc_occ, sizeof(int)*nactive);
    psio_read_entry(CC_INFO, "QT->CC Active Virt Order",
            (char *) moinfo.cc_vir, sizeof(int)*nactive);

  }

  /* Build sosym array (for AO-basis BT2) */
  moinfo.sosym = init_int_array(moinfo.nso);
  for(h=0,q=0; h < nirreps; h++)
    for(p=0; p < moinfo.sopi[h]; p++)
      moinfo.sosym[q++] = h;

  /* Get the active virtual orbitals */
  if(params.ref == 0 || params.ref == 1) { /** RHF/ROHF **/

    C = (double ***) malloc(nirreps * sizeof(double **));
    next = PSIO_ZERO;
    for(h=0; h < nirreps; h++) {
      if(moinfo.sopi[h] && moinfo.virtpi[h]) {
        C[h] = block_matrix(moinfo.sopi[h],moinfo.virtpi[h]);
        psio_read(CC_INFO, "RHF/ROHF Active Virtual Orbitals", (char *) C[h][0],
                  moinfo.sopi[h]*moinfo.virtpi[h]*sizeof(double), next, &next);
      }
    }
    moinfo.C = C;
  }
  else if(params.ref == 2) { /** UHF **/

    Ca = (double ***) malloc(nirreps * sizeof(double **));
    next = PSIO_ZERO;
    for(h=0; h < nirreps; h++) {
      if(moinfo.sopi[h] && moinfo.avirtpi[h]) {
        Ca[h] = block_matrix(moinfo.sopi[h],moinfo.avirtpi[h]);
        psio_read(CC_INFO, "UHF Active Alpha Virtual Orbs", (char *) Ca[h][0],
                  moinfo.sopi[h]*moinfo.avirtpi[h]*sizeof(double), next, &next);
      }
    }
    moinfo.Ca = Ca;


    Cb = (double ***) malloc(nirreps * sizeof(double **));
    next = PSIO_ZERO;
    for(h=0; h < nirreps; h++) {
      if(moinfo.sopi[h] && moinfo.bvirtpi[h]) {
        Cb[h] = block_matrix(moinfo.sopi[h],moinfo.bvirtpi[h]);
        psio_read(CC_INFO, "UHF Active Beta Virtual Orbs", (char *) Cb[h][0],
                  moinfo.sopi[h]*moinfo.bvirtpi[h]*sizeof(double), next, &next);
      }
    }
    moinfo.Cb = Cb;
  }

  /* Compute spatial-orbital reordering arrays */
  if(params.ref == 0 || params.ref == 1) {
    moinfo.pitzer2qt = init_int_array(moinfo.nmo);
    moinfo.qt2pitzer = init_int_array(moinfo.nmo);
    reorder_qt(moinfo.clsdpi, moinfo.openpi, moinfo.frdocc, moinfo.fruocc,
           moinfo.pitzer2qt, moinfo.orbspi, moinfo.nirreps);
    for(i=0; i < moinfo.nmo; i++) {
      j = moinfo.pitzer2qt[i];
      moinfo.qt2pitzer[j] = i;
    }
  }
  else if(params.ref == 2) {
    moinfo.pitzer2qt_a = init_int_array(moinfo.nmo);
    moinfo.qt2pitzer_a = init_int_array(moinfo.nmo);
    moinfo.pitzer2qt_b = init_int_array(moinfo.nmo);
    moinfo.qt2pitzer_b = init_int_array(moinfo.nmo);
    reorder_qt_uhf(moinfo.clsdpi, moinfo.openpi, moinfo.frdocc, moinfo.fruocc,
           moinfo.pitzer2qt_a, moinfo.pitzer2qt_b, moinfo.orbspi,
           moinfo.nirreps);
    for(i=0; i < moinfo.nmo; i++) {
      j = moinfo.pitzer2qt_a[i];
      moinfo.qt2pitzer_a[j] = i;
      j = moinfo.pitzer2qt_b[i];
      moinfo.qt2pitzer_b[j] = i;
    }
  }

  /* Adjust clsdpi array for frozen orbitals */
  for(i=0; i < nirreps; i++)
    moinfo.clsdpi[i] -= moinfo.frdocc[i];

  moinfo.uoccpi = init_int_array(moinfo.nirreps);
  for(i=0; i < nirreps; i++)
    moinfo.uoccpi[i] = moinfo.orbspi[i] - moinfo.clsdpi[i] -
      moinfo.openpi[i] - moinfo.fruocc[i] -
      moinfo.frdocc[i];

  if(params.ref == 0) {
    moinfo.nvirt = 0;
    for(h=0; h < nirreps; h++) moinfo.nvirt += moinfo.virtpi[h];
  }

  psio_read_entry(CC_INFO, "Reference Energy", (char *) &(moinfo.eref),
          sizeof(double));

  fprintf(outfile,"\n\tNuclear Rep. energy (chkpt)   = %20.15f\n",moinfo.enuc);
  fprintf(outfile,  "\tSCF energy          (chkpt)   = %20.15f\n",moinfo.escf);
  fprintf(outfile,  "\tReference energy    (file100) = %20.15f\n",moinfo.eref);
}

/* Frees memory allocated in get_moinfo() and dumps out the energy. */
void cleanup(void)
{
  int i, h;
  char *keyw=NULL;

  /* Save the energy to PSIF_CHKPT as well */
  chkpt_init(PSIO_OPEN_OLD);
  if( params.wfn == "CC2" || params.wfn == "EOM_CC2" ) {
    psio_write_entry(CC_INFO, "CC2 Energy", (char *) &(moinfo.ecc),
             sizeof(double));

    keyw = chkpt_build_keyword("CC2 Energy");
    psio_write_entry(PSIF_CHKPT, keyw, (char *) &(moinfo.ecc),
             sizeof(double));
    free(keyw);
  }
  else if( params.wfn == "CC3"  || params.wfn == "EOM_CC3" ) {
    psio_write_entry(CC_INFO, "CC3 Energy", (char *) &(moinfo.ecc),
             sizeof(double));

    keyw = chkpt_build_keyword("CC3 Energy");
    psio_write_entry(PSIF_CHKPT, keyw, (char *) &(moinfo.ecc),
             sizeof(double));
    free(keyw);
  }
  else {
    psio_write_entry(CC_INFO, "CCSD Energy", (char *) &(moinfo.ecc),
             sizeof(double));

    keyw = chkpt_build_keyword("CCSD Energy");
    psio_write_entry(PSIF_CHKPT, keyw, (char *) &(moinfo.ecc),
             sizeof(double));
    free(keyw);
  }
  chkpt_close();

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
  // Wavefunction owns these arrays
//  free(moinfo.fruocc);
//  free(moinfo.frdocc);
  for(i=0; i < moinfo.nirreps; i++)
    free(moinfo.labels[i]);
  free(moinfo.labels);
  if(params.ref == 2) {
    free(moinfo.aoccpi);
    free(moinfo.boccpi);
    free(moinfo.avirtpi);
    free(moinfo.bvirtpi);
    free(moinfo.aocc_sym);
    free(moinfo.bocc_sym);
    free(moinfo.avir_sym);
    free(moinfo.bvir_sym);
    free(moinfo.aocc_off);
    free(moinfo.bocc_off);
    free(moinfo.avir_off);
    free(moinfo.bvir_off);
    free(moinfo.qt_aocc);
    free(moinfo.qt_bocc);
    free(moinfo.qt_avir);
    free(moinfo.qt_bvir);
    free(moinfo.cc_aocc);
    free(moinfo.cc_bocc);
    free(moinfo.cc_avir);
    free(moinfo.cc_bvir);
  }
  else {
    free(moinfo.occpi);
    free(moinfo.virtpi);
    free(moinfo.occ_sym);
    free(moinfo.vir_sym);
    free(moinfo.occ_off);
    free(moinfo.vir_off);
    free(moinfo.qt_occ);
    free(moinfo.qt_vir);
    free(moinfo.cc_occ);
    free(moinfo.cc_vir);
  }

}

}} // namespace psi::ccenergy
