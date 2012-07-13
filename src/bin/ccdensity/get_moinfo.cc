/*! \file
    \ingroup CCDENSITY
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include <string.h>
#include <string.h>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>
#include <libqt/qt.h>
#include <libchkpt/chkpt.h>
#include "MOInfo.h"
#include "Params.h"
#include "Frozen.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccdensity {

/*
** get_moinfo():  Routine to obtain basic orbital information from
** CHKPT and CC_INFO.
**
** T. Daniel Crawford, October 1996.
** Modified by TDC, March 1999.
*/

void get_moinfo(void)
{
  int i, j, h, errcod;
  int nactive;
  double **scf_pitzer;

  chkpt_init(PSIO_OPEN_OLD);
  moinfo.nirreps = chkpt_rd_nirreps();
  moinfo.nmo = chkpt_rd_nmo();
  moinfo.nso = chkpt_rd_nso();
  moinfo.iopen = chkpt_rd_iopen();
  moinfo.labels = chkpt_rd_irr_labs();
  moinfo.enuc = chkpt_rd_enuc();
  moinfo.escf = chkpt_rd_escf();
  moinfo.orbspi = chkpt_rd_orbspi();
  moinfo.clsdpi = chkpt_rd_clsdpi();
  moinfo.openpi = chkpt_rd_openpi();
  scf_pitzer = chkpt_rd_scf();
  chkpt_close();

  moinfo.sym = 0;
  for (i=0;i<moinfo.nirreps;++i)
    for (j=0;j<moinfo.openpi[i];++j)
      moinfo.sym = moinfo.sym ^ i;

  /* Get frozen and active orbital lookups from CC_INFO */
  moinfo.frdocc = init_int_array(moinfo.nirreps);
  moinfo.fruocc = init_int_array(moinfo.nirreps);
  psio_read_entry(CC_INFO, "Frozen Core Orbs Per Irrep",
		  (char *) moinfo.frdocc, sizeof(int)*moinfo.nirreps);
  psio_read_entry(CC_INFO, "Frozen Virt Orbs Per Irrep",
		  (char *) moinfo.fruocc, sizeof(int)*moinfo.nirreps);

  psio_read_entry(CC_INFO, "No. of Active Orbitals", (char *) &(nactive),
		  sizeof(int)); 
    moinfo.nactive = nactive;

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

  /* Compute spatial-orbital reordering arrays */
  moinfo.pitzer2qt = init_int_array(moinfo.nmo);
  moinfo.qt2pitzer = init_int_array(moinfo.nmo);
  reorder_qt(moinfo.clsdpi, moinfo.openpi, moinfo.frdocc, moinfo.fruocc,
             moinfo.pitzer2qt, moinfo.orbspi, moinfo.nirreps);
  for(i=0; i < moinfo.nmo; i++) {
    j = moinfo.pitzer2qt[i];
    moinfo.qt2pitzer[j] = i;
  }

  /* Adjust clsdpi array for frozen orbitals */
  for(i=0; i < moinfo.nirreps; i++)
    moinfo.clsdpi[i] -= moinfo.frdocc[i];

  moinfo.uoccpi = init_int_array(moinfo.nirreps);
  for(i=0; i < moinfo.nirreps; i++)
    moinfo.uoccpi[i] = moinfo.orbspi[i] - moinfo.clsdpi[i] -
      moinfo.openpi[i] - moinfo.fruocc[i] -
      moinfo.frdocc[i];

  moinfo.nfzc = moinfo.nfzv = moinfo.nclsd = moinfo.nopen = moinfo.nuocc = 0;
  for(h=0; h < moinfo.nirreps; h++) {
    moinfo.nfzc += moinfo.frdocc[h];
    moinfo.nfzv += moinfo.fruocc[h];
    moinfo.nclsd += moinfo.clsdpi[h];
    moinfo.nopen += moinfo.openpi[h];
    moinfo.nuocc += moinfo.uoccpi[h];
  }

  if(params.ref == 0 || params.ref == 1) { /** RHF or ROHF **/

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

    /* Sort SCF MOs to QT order */
    moinfo.scf_qt = block_matrix(moinfo.nmo, moinfo.nmo);
    int I;
    for(i=0; i < moinfo.nmo; i++) {
      I = moinfo.pitzer2qt[i];  /* Pitzer --> QT */
      for(j=0; j < moinfo.nmo; j++) moinfo.scf_qt[j][I] = scf_pitzer[j][i];
    }
    free_block(scf_pitzer);
  }
  else if(params.ref == 2) { /** UHF **/

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

  psio_read_entry(CC_INFO, "Reference Energy", (char *) &(moinfo.eref),
		  sizeof(double));

  fprintf(outfile,"\n\tNuclear Rep. energy (chkpt)   = %20.15f\n",moinfo.enuc);
  fprintf(outfile,  "\tSCF energy          (chkpt)   = %20.15f\n",moinfo.escf);
  fprintf(outfile,  "\tReference energy    (file100) = %20.15f\n",moinfo.eref);

  if(params.wfn == "CC2" || params.wfn == "EOM_CC2") {
    psio_read_entry(CC_INFO, "CC2 Energy", (char *) &(moinfo.ecc),
                    sizeof(double));
    fprintf(outfile,  "\tCC2 energy          (CC_INFO) = %20.15f\n",moinfo.ecc);
    fprintf(outfile,  "\tTotal CC2 energy    (CC_INFO) = %20.15f\n",
            moinfo.eref+moinfo.ecc);
  }
  else if( params.wfn == "CCSD" || params.wfn == "EOM_CCSD") {
    psio_read_entry(CC_INFO, "CCSD Energy", (char *) &(moinfo.ecc),
                    sizeof(double));
    fprintf(outfile,  "\tCCSD energy         (CC_INFO) = %20.15f\n",moinfo.ecc);
    fprintf(outfile,  "\tTotal CCSD energy   (CC_INFO) = %20.15f\n",
            moinfo.eref+moinfo.ecc);
  }
  else if(params.wfn == "CCSD_T") {
    psio_read_entry(CC_INFO, "CCSD Energy", (char *) &(moinfo.ecc), sizeof(double));
    psio_read_entry(CC_INFO, "(T) Energy", (char *) &(moinfo.et), sizeof(double));
    fprintf(outfile,  "\tCCSD energy         (CC_INFO) = %20.15f\n",moinfo.ecc);
    fprintf(outfile,  "\t(T) energy          (CC_INFO) = %20.15f\n",moinfo.et);
    fprintf(outfile,  "\tTotal CCSD(T) energy(CC_INFO) = %20.15f\n",
            moinfo.eref+moinfo.ecc+moinfo.et);
  }
  else if(params.wfn == "CC3" || params.wfn == "EOM_CC3") {
    psio_read_entry(CC_INFO, "CC3 Energy", (char *) &(moinfo.ecc),
                    sizeof(double));
    fprintf(outfile,  "\tCC3 energy          (CC_INFO) = %20.15f\n",moinfo.ecc);
    fprintf(outfile,  "\tTotal CC3 energy    (CC_INFO) = %20.15f\n",
            moinfo.eref+moinfo.ecc);
  }

  fflush(outfile);
}

/* Frees memory allocated in get_moinfo(). */
void cleanup(void)
{
  int i;

  free(moinfo.orbspi);
  free(moinfo.clsdpi);
  free(moinfo.openpi);
  free(moinfo.uoccpi);
  free(moinfo.fruocc);
  free(moinfo.frdocc);
  for(i=0; i < moinfo.nirreps; i++)
    free(moinfo.labels[i]);
  free(moinfo.labels);

  if(params.ref == 0 || params.ref == 1) { /** RHF or ROHF **/
    free(moinfo.occ_sym);
    free(moinfo.vir_sym);
    free(moinfo.occ_off);
    free(moinfo.vir_off);
    free(moinfo.occpi);
    free(moinfo.virtpi);
    free(moinfo.qt_occ);
    free(moinfo.qt_vir);
    free(moinfo.cc_occ);
    free(moinfo.cc_vir);
    free(moinfo.pitzer2qt);
    free(moinfo.qt2pitzer);
    free_block(moinfo.scf_qt);
  }
  else if(params.ref == 2) { /** UHF **/
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
    free(moinfo.qt_aocc);
    free(moinfo.qt_bocc);
    free(moinfo.qt_avir);
    free(moinfo.qt_bvir);
    free(moinfo.cc_aocc);
    free(moinfo.cc_bocc);
    free(moinfo.cc_avir);
    free(moinfo.cc_bvir);
  }
}


}} // namespace psi::ccdensity
