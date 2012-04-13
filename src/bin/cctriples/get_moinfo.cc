/*! \file
    \ingroup CCTRIPLES
    \brief Enter brief description of file here
*/
#include <cstdio>
#include <cstdlib>
#include <string>
#include <cstring>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>
#include <libchkpt/chkpt.h>
#include <liboptions/liboptions.h>
#include "psi4-dec.h"
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace cctriples {

    /*
    ** get_moinfo():  Routine to obtain basic orbital information from
    ** CHKPT and CC_INFO.
    **
    ** T. Daniel Crawford, October 1996.
    ** Modified by TDC, March 1999.
    */

    void get_moinfo(Options &options)
    {
      int i, h, errcod, nactive, nirreps;
      std::string junk;
      chkpt_init(PSIO_OPEN_OLD);
      moinfo.nirreps = chkpt_rd_nirreps();
      moinfo.nmo = chkpt_rd_nmo();
      moinfo.iopen = chkpt_rd_iopen();
      moinfo.labels = chkpt_rd_irr_labs();
      moinfo.enuc = chkpt_rd_enuc();
      moinfo.escf = chkpt_rd_escf();
      moinfo.orbspi = chkpt_rd_orbspi();
      moinfo.clsdpi = chkpt_rd_clsdpi();
      moinfo.openpi = chkpt_rd_openpi();
      moinfo.phase = chkpt_rd_phase_check();
      chkpt_close();

      nirreps = moinfo.nirreps;

      params.wfn = options.get_str("WFN");
      if( params.wfn != "CCSD" && params.wfn != "CCSD_T" &&
                params.wfn != "BCCD" && params.wfn != "BCCD_T" ) {
        throw PsiException("Invalid value of input keyword WFN",__FILE__,__LINE__);
      }

  params.nthreads = Process::environment.get_n_threads();
  if (options["CC_NUM_THREADS"].has_changed()){
     params.nthreads = options.get_int("CC_NUM_THREADS");
  }

  params.semicanonical = 0;
  junk = options.get_str("REFERENCE");
  /* if no reference is given, assume rhf */
  if(junk == "RHF") params.ref = 0;
  else if(junk == "ROHF" && params.wfn == "CCSD_T") {
    params.ref = 2;
    params.semicanonical = 1;
  }
  else if(junk == "ROHF") params.ref = 1;
  else if(junk == "UHF") params.ref = 2;
  else {
    throw PsiException("Invalid value of input keyword REFERENCE",__FILE__,__LINE__);
  }

  junk = options.get_str("DERTYPE");
  if(junk == "NONE") params.dertype = 0;
  else if(junk == "FIRST") params.dertype = 1;
  else {
          throw PsiException("Value of keyword DERTYPE is not applicable to CCSD(T)",__FILE__,__LINE__);
        }

      /* Get frozen and active orbital lookups from CC_INFO */
      moinfo.frdocc = init_int_array(moinfo.nirreps);
      moinfo.fruocc = init_int_array(moinfo.nirreps);
      psio_read_entry(CC_INFO, "Frozen Core Orbs Per Irrep",
                      (char *) moinfo.frdocc, sizeof(int)*moinfo.nirreps);
      psio_read_entry(CC_INFO, "Frozen Virt Orbs Per Irrep",
                      (char *) moinfo.fruocc, sizeof(int)*moinfo.nirreps);

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

      }
      else { /** RHF or ROHF **/

        moinfo.occpi = init_int_array(moinfo.nirreps);
        moinfo.virtpi = init_int_array(moinfo.nirreps);
        psio_read_entry(CC_INFO, "Active Occ Orbs Per Irrep",
                        (char *) moinfo.occpi, sizeof(int)*moinfo.nirreps);
        psio_read_entry(CC_INFO, "Active Virt Orbs Per Irrep",
                        (char *) moinfo.virtpi, sizeof(int)*moinfo.nirreps);

        psio_read_entry(CC_INFO, "No. of Active Orbitals", (char *) &(nactive),
                        sizeof(int));

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

      /* Adjust clsdpi array for frozen orbitals */
      for(i=0; i < moinfo.nirreps; i++)
        moinfo.clsdpi[i] -= moinfo.frdocc[i];

      moinfo.uoccpi = init_int_array(moinfo.nirreps);
      for(i=0; i < moinfo.nirreps; i++)
        moinfo.uoccpi[i] = moinfo.orbspi[i] - moinfo.clsdpi[i] -
          moinfo.openpi[i] - moinfo.fruocc[i] -
          moinfo.frdocc[i];

      fprintf(outfile,"\n\n");
      fprintf(outfile, "\tWave function   =    %6s\n",params.wfn.c_str());
      if(params.semicanonical) {
        fprintf(outfile, "\tReference wfn   =    ROHF changed to UHF for Semicanonical Orbitals\n");
      }
      else {
        fprintf(outfile, "\tReference wfn   =    %5s\n",
                (params.ref == 0) ? "RHF" : ((params.ref == 1) ? "ROHF" : "UHF"));
      }
      psio_read_entry(CC_INFO, "Reference Energy", (char *) &(moinfo.eref),
                      sizeof(double));
      psio_read_entry(CC_INFO, "CCSD Energy", (char *) &(moinfo.ecc),
                      sizeof(double));

      fprintf(outfile,"\n\tNuclear Rep. energy (chkpt)   = %20.15f\n",moinfo.enuc);
      fprintf(outfile,  "\tSCF energy          (chkpt)   = %20.15f\n",moinfo.escf);
      fprintf(outfile,  "\tReference energy    (file100) = %20.15f\n",moinfo.eref);
      fprintf(outfile,  "\tCCSD energy         (file100) = %20.15f\n",moinfo.ecc);
      fprintf(outfile,  "\tTotal CCSD energy   (file100) = %20.15f\n",
              moinfo.eref+moinfo.ecc);
    }

    /* Frees memory allocated in get_moinfo() and dumps some info. */
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
      if(params.ref == 2) {
        free(moinfo.aoccpi);
        free(moinfo.boccpi);
        free(moinfo.avirtpi);
        free(moinfo.bvirtpi);
        free(moinfo.aocc_sym);
        free(moinfo.bocc_sym);
        free(moinfo.avir_sym);
        free(moinfo.bvir_sym);
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


  }} // namespace psi::CCTRIPLES
