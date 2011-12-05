/*! \file
    \ingroup CCTRIPLES
    \brief Enter brief description of file here
*/
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>
#include <libchkpt/chkpt.h>
#include <libdpd/dpd.h>
#include <libqt/qt.h>
#include <psifiles.h>
#include <psi4-dec.h>
#include "Params.h"
#include "MOInfo.h"
#include "globals.h"

namespace psi { namespace cctriples {


    void init_io();
    void title(void);
    void get_moinfo(Options&);
    void exit_io(void);
    void cleanup(void);
    double ET_RHF(void);
    double ET_AAA(void);
    double ET_AAB(void);
    double ET_ABB(void);
    double ET_BBB(void);
    double ET_UHF_AAA(void);
    double ET_UHF_BBB(void);
    double ET_UHF_AAB(void);
    double ET_UHF_ABB(void);
    void count_ijk(void);
    void setup(void);
    int **cacheprep_rhf(int level, int *cachefiles);
    int **cacheprep_uhf(int level, int *cachefiles);
    void cachedone_uhf(int **cachelist);
    void cachedone_rhf(int **cachelist);

    void T3_grad_RHF(void);
    void T3_grad_UHF_AAA(void);
    void T3_grad_UHF_BBB(void);
    void T3_grad_UHF_AAB(void);
    void T3_grad_UHF_BBA(void);


    void T3_UHF_AAA_abc(double ***W, double ***V, int disc, int nirreps, int A, int Ga, int B, int Gb, int C, int Gc,
            dpdbuf4 *C2, dpdbuf4 *F, dpdbuf4 *E, dpdfile2 *C1, dpdbuf4 *D, dpdfile2 *fIA, dpdfile2 *fIJ, dpdfile2 *fAB,
            int *occpi, int *occ_off, int *virtpi, int *vir_off, double omega);

    void T3_UHF_AAB_abc(double ***W, double ***V, int disc, int nirreps,
            int I, int Gi, int J, int Gj, int K, int Gk,
            dpdbuf4 *T2AA, dpdbuf4 *T2AB, dpdbuf4 *T2BA, dpdbuf4 *FAA, dpdbuf4 *FAB, dpdbuf4 *FBA,
            dpdbuf4 *EAA, dpdbuf4 *EAB, dpdbuf4 *EBA, dpdfile2 *T1A, dpdfile2 *T1B,
            dpdbuf4 *DAA, dpdbuf4 *DAB, dpdfile2 *fIA, dpdfile2 *fia,
            dpdfile2 *fIJ, dpdfile2 *fij,dpdfile2 *fAB, dpdfile2 *fab,
            int *aoccpi, int *aocc_off, int *boccpi, int *bocc_off,
            int *avirtpi, int *avir_off, int *bvirtpi, int *bvir_off, double omega);
    void transpose_integrals();
    void test_abc_loops_AAA();
    void test_abc_loops_AAB();
    void test_abc_loops_BBA();
    void test_abc_loops_BBB();

PsiReturnType cctriples(Options &options)
{
  double ETAAA, ETAAB, ETABB, ETBBB, ET;
  long int memory;
  int **cachelist, *cachefiles;
  dpdfile2 T1;
  double **geom, *zvals, value;
  FILE *efile;
  int i, errcod, natom;
  char *keyw = NULL;

  init_io();
  title();

  timer_on("CCtriples");

  get_moinfo(options);
  memory = Process::environment.get_memory();

  cachefiles = init_int_array(PSIO_MAXUNIT);


  if(params.ref == 0) { /*** RHF ***/
    cachelist = cacheprep_rhf(2, cachefiles);

    dpd_init(0, moinfo.nirreps, memory, 0, cachefiles, cachelist, NULL,
         2, moinfo.occpi, moinfo.occ_sym, moinfo.virtpi, moinfo.vir_sym);
  }
  else if(params.ref == 2) { /*** UHF ***/
    cachelist = cacheprep_uhf(2, cachefiles);

    dpd_init(0, moinfo.nirreps, memory, 0, cachefiles,
         cachelist, NULL, 4, moinfo.aoccpi, moinfo.aocc_sym, moinfo.avirtpi,
         moinfo.avir_sym, moinfo.boccpi, moinfo.bocc_sym, moinfo.bvirtpi, moinfo.bvir_sym);
  }

  count_ijk();
  fflush(outfile);

  if(params.ref == 0) { /** RHF **/

    ET = ET_RHF();
    fprintf(outfile, "\t(T) energy                    = %20.15f\n", ET);
    fprintf(outfile, "      * CCSD(T) total energy          = %20.15f\n",
        ET + moinfo.ecc + moinfo.eref);

    Process::environment.globals["(T) CORRECTION ENERGY"] = ET;
    Process::environment.globals["CCSD(T) CORRELATION ENERGY"] = ET + moinfo.ecc;
    Process::environment.globals["CCSD(T) TOTAL ENERGY"] = ET + moinfo.ecc + moinfo.eref;

    /* Compute triples contributions to the gradient */
    if(params.dertype == 1){
      T3_grad_RHF();
    }
  }
  else if(params.ref == 1 ) { /** ROHF --- don't use this right now! **/

    throw PsiException("ROHF-CCSD(T) is not yet available",__FILE__,__LINE__);

    ETAAA = ET_AAA();
    fprintf(outfile, "\tAAA (T) energy                = %20.15f\n", ETAAA);
    ETAAB = ET_AAB();
    fprintf(outfile, "\tAAB (T) energy                = %20.15f\n", ETAAB);
    ETABB = ET_ABB();
    fprintf(outfile, "\tABB (T) energy                = %20.15f\n", ETABB);
    ETBBB = ET_BBB();
    fprintf(outfile, "\tBBB (T) energy                = %20.15f\n", ETBBB);
    ET = ETAAA + ETAAB + ETABB + ETBBB;
    fprintf(outfile, "\t(T) energy                    = %20.15f\n", ET);
    fprintf(outfile, "      * CCSD(T) total energy          = %20.15f\n",
        ET + moinfo.ecc + moinfo.eref);

    Process::environment.globals["AAA (T) CORRECTION ENERGY"] = ETAAA;
    Process::environment.globals["AAB (T) CORRECTION ENERGY"] = ETAAB;
    Process::environment.globals["ABB (T) CORRECTION ENERGY"] = ETABB;
    Process::environment.globals["BBB (T) CORRECTION ENERGY"] = ETBBB;
    Process::environment.globals["(T) CORRECTION ENERGY"] = ET;
    Process::environment.globals["CCSD(T) CORRELATION ENERGY"] = ET + moinfo.ecc;
    Process::environment.globals["CCSD(T) TOTAL ENERGY"] = ET + moinfo.ecc + moinfo.eref;
  }
  else if(params.ref == 2) { /** UHF **/
    ETAAA = ET_UHF_AAA();
    fprintf(outfile, "\tAAA (T) energy                = %20.15f\n", ETAAA);
    fflush(outfile);

    ETBBB = ET_UHF_BBB();
    fprintf(outfile, "\tBBB (T) energy                = %20.15f\n", ETBBB);
    fflush(outfile);

    ETAAB = ET_UHF_AAB();
    fprintf(outfile, "\tAAB (T) energy                = %20.15f\n", ETAAB);
    fflush(outfile);

    ETABB = ET_UHF_ABB();
    fprintf(outfile, "\tABB (T) energy                = %20.15f\n", ETABB);
    fflush(outfile);

    ET = ETAAA + ETAAB + ETABB + ETBBB;
    fprintf(outfile, "\t(T) energy                    = %20.15f\n", ET);
    fprintf(outfile, "      * CCSD(T) total energy          = %20.15f\n",
        ET + moinfo.ecc + moinfo.eref);

    Process::environment.globals["AAA (T) CORRECTION ENERGY"] = ETAAA;
    Process::environment.globals["AAB (T) CORRECTION ENERGY"] = ETAAB;
    Process::environment.globals["ABB (T) CORRECTION ENERGY"] = ETABB;
    Process::environment.globals["BBB (T) CORRECTION ENERGY"] = ETBBB;
    Process::environment.globals["(T) CORRECTION ENERGY"] = ET;
    Process::environment.globals["CCSD(T) CORRELATION ENERGY"] = ET + moinfo.ecc;
    Process::environment.globals["CCSD(T) TOTAL ENERGY"] = ET + moinfo.ecc + moinfo.eref;

    if(params.dertype==1) {

      transpose_integrals();
      test_abc_loops_AAA();
      test_abc_loops_BBB();
      test_abc_loops_AAB();
      test_abc_loops_BBA();

      fprintf(outfile, "\n\tComputing (T) contributions to CC density...\n");
      fflush(outfile);
      T3_grad_UHF_AAA();
      fprintf(outfile, "\tAAA contributions complete.\n");
      fflush(outfile);
      T3_grad_UHF_BBB();
      fprintf(outfile, "\tBBB contributions complete.\n");
      fflush(outfile);
      T3_grad_UHF_AAB();
      fprintf(outfile, "\tAAB contributions complete.\n");
      fflush(outfile);
      T3_grad_UHF_BBA();
      fprintf(outfile, "\tBBA contributions complete.\n");
      fflush(outfile);
    }
  }

  fprintf(outfile, "\n");

  /* Write total energy and (T) contribution to the checkpoint file */
  chkpt_init(PSIO_OPEN_OLD);
  chkpt_wt_etot(ET+moinfo.ecc+moinfo.eref);
  chkpt_wt_e_t(ET);
  chkpt_close();

  /* Write pertinent data to energy.dat */
//  if(params.wfn == "CCSD_T" || params.wfn == "BCCD_T") {
//    chkpt_init(PSIO_OPEN_OLD);
//    natom = chkpt_rd_natom();
//    geom = chkpt_rd_geom();
//    zvals = chkpt_rd_zvals();
//    chkpt_close();
//    ffile(&efile,"energy.dat",1);
//    fprintf(efile, "*\n");
//    for(i=0; i < natom; i++)
//      fprintf(efile, " %4d   %5.2f     %13.10f    %13.10f    %13.10f\n",
//          i+1, zvals[i], geom[i][0], geom[i][1], geom[i][2]);
//    free_block(geom);  free(zvals);
//    fprintf(efile, "SCF(30)   %22.12f\n", moinfo.escf);
//    fprintf(efile, "REF(100)  %22.12f\n", moinfo.eref);
//    if(params.wfn == "CCSD_T") {
//      fprintf(efile, "CCSD      %22.12f\n", (moinfo.ecc+moinfo.eref));
//      fprintf(efile, "CCSD(T)   %22.12f\n", (ET+ moinfo.ecc+moinfo.eref));
//    }
//    else if(params.wfn == "BCCD_T") {
//      fprintf(efile, "BCCD      %22.12f\n", (moinfo.ecc+moinfo.eref));
//      fprintf(efile, "BCCD(T)   %22.12f\n", (ET+ moinfo.ecc+moinfo.eref));
//    }
//    fclose(efile);
//  }

  /* Dump triples energy to CC_INFO */
  psio_write_entry(CC_INFO, "(T) Energy", (char *) &(ET), sizeof(double));

  Process::environment.globals["(T) CORRECTION ENERGY"] = ET; 
  Process::environment.globals["CCSD(T) TOTAL ENERGY"] = ET+ moinfo.ecc+moinfo.eref;
  Process::environment.globals["CCSD(T) CORRELATION ENERGY"] = ET+ moinfo.ecc;
  Process::environment.globals["CURRENT ENERGY"] = ET+ moinfo.ecc+moinfo.eref;
  Process::environment.globals["CURRENT CORRELATION ENERGY"] = ET+ moinfo.ecc;

  dpd_close(0);

  if(params.ref == 2) cachedone_uhf(cachelist);
  else cachedone_rhf(cachelist);
  free(cachefiles);

  cleanup();

  timer_off("CCtriples");

  exit_io();
  return Success;
}

void init_io()
{
  tstart();

  for(int i=CC_MIN; i <= CC_MAX; i++) psio_open(i,1);
}

void title(void)
{
  fprintf(outfile, "\t\t\t**************************\n");
  fprintf(outfile, "\t\t\t*                        *\n");
  fprintf(outfile, "\t\t\t*        CCTRIPLES       *\n");
  fprintf(outfile, "\t\t\t*                        *\n");
  fprintf(outfile, "\t\t\t**************************\n");
}

void exit_io(void)
{
  int i;
  for(i=CC_MIN; i <= CC_MAX; i++) psio_close(i,1);
  tstop();
}


}} // namespace psi::CCTRIPLES
