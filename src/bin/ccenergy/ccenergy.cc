/*! \file
    \ingroup CCENERGY
    \brief Enter brief description of file here
*/
/*
**  CCENERGY: Program to calculate coupled cluster energies.
*/

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <sys/types.h>
#include <unistd.h>
#include <libciomr/libciomr.h>
#include <libdpd/dpd.h>
#include <libchkpt/chkpt.h>
#include <libqt/qt.h>
#include <libint/libint.h>
#include <libmints/wavefunction.h>
#include <libpsio/psio.h>
#include <libpsio/psio.hpp>
#include <sys/types.h>
#include <psifiles.h>
#include "Params.h"
#include "MOInfo.h"
#include "Local.h"
#include "globals.h"
#include "ccwave.h"

namespace psi { namespace ccenergy {

#define IOFF_MAX 32641

void init_io();
void init_ioff(void);
void title(void);
void get_moinfo(void);
void get_params(Options &);
void init_amps(void);
void tau_build(void);
void taut_build(void);
double energy(void);
double mp2_energy(void);
void sort_amps(void);
void Fae_build(void);
void Fmi_build(void);
void Fme_build(void);
void t1_build(void);
void Wmnij_build(void);
void Z_build(void);
void Y_build(void);
void X_build(void);
void Wmbej_build(void);
void t2_build(void);
void tsave(void);
int converged(void);
double diagnostic(void);
double d1diag(void);
double new_d1diag(void);
void exit_io(void);
void cleanup(void);
void update(void);
void diis(int iter);
void ccdump(void);
int **cacheprep_uhf(int level, int *cachefiles);
int **cacheprep_rhf(int level, int *cachefiles);
void cachedone_rhf(int **cachelist);
void cachedone_uhf(int **cachelist);
void memchk(void);
struct dpd_file4_cache_entry *priority_list(void);
void spinad_amps(void);
void status(const char *, FILE *);
void lmp2(void);
void amp_write(void);
void amp_write(void);
int rotate(void);
double **fock_build(double **D);
void analyze(void);
void cc3_Wmnie(void);
void cc3_Wamef(void);
void cc3_Wmnij(void);
void cc3_Wmbij(void);
void cc3_Wabei(void);
void cc3(void);
void cc2_Wmnij_build(void);
void cc2_Wmbij_build(void);
void cc2_Wabei_build(void);
void cc2_t2_build(void);
void one_step(void);
void denom(void);
void pair_energies(double** epair_aa, double** epair_ab);
void print_pair_energies(double* emp2_aa, double* emp2_ab, double* ecc_aa,
double* ecc_ab);
void checkpoint(void);

/* local correlation functions */
void local_init(void);
void local_done(void);

PsiReturnType ccenergy(Options &options);

}} //namespace psi::ccenergy

// Forward declaration to call cctriples
namespace psi { namespace cctriples {
PsiReturnType cctriples(Options &options);
}}

namespace psi { namespace ccenergy {

CCEnergyWavefunction::CCEnergyWavefunction(boost::shared_ptr<Wavefunction> reference_wavefunction, Options &options)
    : Wavefunction(options, _default_psio_lib_)
{
    set_reference_wavefunction(reference_wavefunction);
    init();
}

CCEnergyWavefunction::~CCEnergyWavefunction()
{

}

void CCEnergyWavefunction::init()
{
    // Wavefunction creates a chkpt object for you, but we're not going to use it.
    // Destroy it. Otherwise we will see a "file already open" error.
    chkpt_.reset();

    nso_        = reference_wavefunction_->nso();
    nirrep_     = reference_wavefunction_->nirrep();
    nmo_        = reference_wavefunction_->nmo();
    for(int h = 0; h < nirrep_; ++h){
        soccpi_[h] = reference_wavefunction_->soccpi()[h];
        doccpi_[h] = reference_wavefunction_->doccpi()[h];
        frzcpi_[h] = reference_wavefunction_->frzcpi()[h];
        frzvpi_[h] = reference_wavefunction_->frzvpi()[h];
        nmopi_[h]  = reference_wavefunction_->nmopi()[h];
        nsopi_[h]  = reference_wavefunction_->nsopi()[h];
    }
}

double CCEnergyWavefunction::compute_energy()
{
    energy_ = 0.0;
    PsiReturnType ccsd_return;
    if ((ccsd_return = psi::ccenergy::ccenergy(options_)) == Success) {
        // Get the total energy of the CCSD wavefunction
        energy_ = Process::environment.globals["CURRENT ENERGY"];
    }

    if (options_.get_str("WFN") == "CCSD_T") {
        // Make sure ccenergy returned Success
        if (ccsd_return != Success)
            throw PSIEXCEPTION("CCEnergyWavefunction: CCSD did not converge, will not proceed to (T) correction.");

        // Run cctriples
        if (psi::cctriples::cctriples(options_) == Success)
            energy_ = Process::environment.globals["CURRENT ENERGY"];
        else
            energy_ = 0.0;
    }

    return energy_;
}

PsiReturnType ccenergy(Options &options)
{
  int done=0, brueckner_done=0;
  int h, i, j, a, b, row, col, natom;
  double **geom, *zvals, value;
  FILE *efile;
  int **cachelist, *cachefiles;
  struct dpd_file4_cache_entry *priority;
  dpdfile2 t1;
  dpdbuf4 t2;
  double *emp2_aa, *emp2_ab, *ecc_aa, *ecc_ab, tval;

  moinfo.iter=0;

  init_io();
  init_ioff();
  title();

#ifdef TIME_CCENERGY
  timer_on("CCEnergy");
#endif

  get_moinfo();
  get_params(options);

  cachefiles = init_int_array(PSIO_MAXUNIT);

  if(params.ref == 2) { /** UHF **/
    cachelist = cacheprep_uhf(params.cachelev, cachefiles);

    dpd_init(0, moinfo.nirreps, params.memory, 0, cachefiles,
         cachelist, NULL, 4, moinfo.aoccpi, moinfo.aocc_sym, moinfo.avirtpi,
         moinfo.avir_sym, moinfo.boccpi, moinfo.bocc_sym, moinfo.bvirtpi, moinfo.bvir_sym);

    if( params.aobasis != "NONE" ) { /* Set up new DPD's for AO-basis algorithm */
      dpd_init(1, moinfo.nirreps, params.memory, 0, cachefiles, cachelist, NULL,
               4, moinfo.aoccpi, moinfo.aocc_sym, moinfo.sopi, moinfo.sosym,
               moinfo.boccpi, moinfo.bocc_sym, moinfo.sopi, moinfo.sosym);
      dpd_set_default(0);
    }

  }
  else { /** RHF or ROHF **/
    cachelist = cacheprep_rhf(params.cachelev, cachefiles);

    priority = priority_list();

    dpd_init(0, moinfo.nirreps, params.memory, params.cachetype, cachefiles,
         cachelist, priority, 2, moinfo.occpi, moinfo.occ_sym,
         moinfo.virtpi, moinfo.vir_sym);

    if( params.aobasis != "NONE") { /* Set up new DPD for AO-basis algorithm */
      dpd_init(1, moinfo.nirreps, params.memory, 0, cachefiles, cachelist, NULL,
               2, moinfo.occpi, moinfo.occ_sym, moinfo.sopi, moinfo.sosym);
      dpd_set_default(0);
    }

  }

  if ( (params.just_energy) || (params.just_residuals) ) {
    one_step();
    if(params.ref == 2) cachedone_uhf(cachelist); else cachedone_rhf(cachelist);
    free(cachefiles);
    cleanup();
    exit_io();
    return Success;
  }

  if(params.local) {
    local_init();
    if(local.weakp=="MP2") lmp2();
  }

  init_amps();

  /* Compute the MP2 energy while we're here */
  if(params.ref == 0 || params.ref == 2) {
    moinfo.emp2 = mp2_energy();
    psio_write_entry(CC_INFO, "MP2 Energy", (char *) &(moinfo.emp2),sizeof(double));
    Process::environment.globals["MP2 CORRELATION ENERGY"] = moinfo.emp2;
    Process::environment.globals["MP2 TOTAL ENERGY"] = moinfo.emp2 + moinfo.eref;
  }

  if(params.print_mp2_amps) amp_write();

  tau_build();
  taut_build();
  fprintf(outfile, "\t            Solving CC Amplitude Equations\n");
  fprintf(outfile, "\t            ------------------------------\n");
  fprintf(outfile, "  Iter             Energy              RMS        T1Diag      D1Diag    New D1Diag\n");
  fprintf(outfile, "  ----     ---------------------    ---------   ----------  ----------  ----------\n");
  moinfo.ecc = energy();
  pair_energies(&emp2_aa, &emp2_ab);

  moinfo.t1diag = diagnostic();
  moinfo.d1diag = d1diag();
  moinfo.new_d1diag = new_d1diag();
  update();
  checkpoint();
  for(moinfo.iter=1; moinfo.iter <= params.maxiter; moinfo.iter++) {

    sort_amps();

#ifdef TIME_CCENERGY
    timer_on("F build");
#endif
    Fme_build(); Fae_build(); Fmi_build();
    if(params.print & 2) status("F intermediates", outfile);
#ifdef TIME_CCENERGY
    timer_off("F build");
#endif

    t1_build();
    if(params.print & 2) status("T1 amplitudes", outfile);

    if( params.wfn == "CC2"  || params.wfn == "EOM_CC2" ) {

      cc2_Wmnij_build();
      if(params.print & 2) status("Wmnij", outfile);

#ifdef TIME_CCENERGY
      timer_on("Wmbij build");
#endif
      cc2_Wmbij_build();
      if(params.print & 2) status("Wmbij", outfile);
#ifdef TIME_CCENERGY
      timer_off("Wmbij build");
#endif

#ifdef TIME_CCENERGY
      timer_on("Wabei build");
#endif
      cc2_Wabei_build();
      if(params.print & 2) status("Wabei", outfile);
#ifdef TIME_CCENERGY
      timer_off("Wabei build");
#endif

#ifdef TIME_CCENERGY
      timer_on("T2 Build");
#endif
      cc2_t2_build();
      if(params.print & 2) status("T2 amplitudes", outfile);
#ifdef TIME_CCENERGY
      timer_off("T2 Build");
#endif

    }

    else {

#ifdef TIME_CCENERGY
      timer_on("Wmbej build");
#endif
      Wmbej_build();
      if(params.print & 2) status("Wmbej", outfile);
#ifdef TIME_CCENERGY
      timer_off("Wmbej build");
#endif

      Z_build();
      if(params.print & 2) status("Z", outfile);
      Wmnij_build();
      if(params.print & 2) status("Wmnij", outfile);

#ifdef TIME_CCENERGY
      timer_on("T2 Build");
#endif
      t2_build();
      if(params.print & 2) status("T2 amplitudes", outfile);
#ifdef TIME_CCENERGY
      timer_off("T2 Build");
#endif

      if( params.wfn == "CC3" || params.wfn == "EOM_CC3" ) {

        /* step1: build cc3 intermediates, Wabei, Wmnie, Wmbij, Wamef */
        cc3_Wmnij();
        cc3_Wmbij();
        cc3_Wmnie();
        cc3_Wamef();
        cc3_Wabei();

        /* step2: loop over T3's and add contributions to T1 and T2 as you go */
        cc3();
      }
    }

    if (!params.just_residuals)
      denom(); /* apply denominators to T1 and T2 */

    if(converged()) {
      done = 1;

      tsave();
      tau_build(); taut_build();
      moinfo.ecc = energy();
      moinfo.t1diag = diagnostic();
      moinfo.d1diag = d1diag();
      moinfo.new_d1diag = new_d1diag();
      sort_amps();
      update();
      fprintf(outfile, "\n\tIterations converged.\n");
      fflush(outfile);
      fprintf(outfile, "\n");
      amp_write();
      if (params.analyze != 0) analyze();
      break;
    }
    if(params.diis) diis(moinfo.iter);
    tsave();
    tau_build(); taut_build();
    moinfo.ecc = energy();
    moinfo.t1diag = diagnostic();
    moinfo.d1diag = d1diag();
    moinfo.new_d1diag = new_d1diag();
    update();
    checkpoint();
  }
  fprintf(outfile, "\n");
  if(!done) {
    fprintf(outfile, "\t ** Wave function not converged to %2.1e ** \n",
        params.convergence);
    fflush(outfile);
    if( params.aobasis != "NONE" ) dpd_close(1);
    dpd_close(0);
    cleanup();
#ifdef TIME_CCENERGY
    timer_off("CCEnergy");
#endif
    exit_io();
    return Failure;
  }

  fprintf(outfile, "\tSCF energy       (chkpt)              = %20.15f\n", moinfo.escf);
  fprintf(outfile, "\tReference energy (file100)            = %20.15f\n", moinfo.eref);

  //Process::environment.globals["SCF TOTAL ENERGY (CHKPT)"] = moinfo.escf;
  //Process::environment.globals["SCF TOTAL ENERGY"] = moinfo.eref;

  if(params.ref == 0 || params.ref == 2) {
    if (params.scs) {
      fprintf(outfile, "\n\tOS SCS-MP2 correlation energy      = %20.15f\n", moinfo.emp2_os*params.scsmp2_scale_os);
      fprintf(outfile, "\tSS SCS-MP2 correlation energy      = %20.15f\n", moinfo.emp2_ss*params.scsmp2_scale_ss);
      fprintf(outfile, "\tSCS-MP2 correlation energy            = %20.15f\n", moinfo.emp2_os*params.scsmp2_scale_os
            + moinfo.emp2_ss*params.scsmp2_scale_ss);
      fprintf(outfile, "      * SCS-MP2 total energy                  = %20.15f\n", moinfo.eref
            + moinfo.emp2_os*params.scsmp2_scale_os + moinfo.emp2_ss*params.scsmp2_scale_ss);

      Process::environment.globals["SCS-MP2 OPPOSITE-SPIN CORRELATION ENERGY"] = moinfo.emp2_os*params.scsmp2_scale_os;
      Process::environment.globals["SCS-MP2 SAME-SPIN CORRELATION ENERGY"] = moinfo.emp2_ss*params.scsmp2_scale_ss;
      Process::environment.globals["SCS-MP2 CORRELATION ENERGY"] = moinfo.emp2_os*params.scsmp2_scale_os +
            moinfo.emp2_ss*params.scsmp2_scale_ss;
      Process::environment.globals["SCS-MP2 TOTAL ENERGY"] = moinfo.eref +
            moinfo.emp2_os*params.scsmp2_scale_os + moinfo.emp2_ss*params.scsmp2_scale_ss;
    }
    if (params.scsn) {
      fprintf(outfile, "\n\tOS SCSN-MP2 correlation energy      = %20.15f\n", 0.0);
      fprintf(outfile, "\tSS SCSN-MP2 correlation energy      = %20.15f\n", moinfo.emp2_ss*1.76);
      fprintf(outfile, "\tSCSN-MP2 correlation energy            = %20.15f\n", moinfo.emp2_ss*1.76);
      fprintf(outfile, "      * SCSN-MP2 total energy                  = %20.15f\n", moinfo.eref + moinfo.emp2_ss*1.76);

      Process::environment.globals["SCSN-MP2 OPPOSITE-SPIN CORRELATION ENERGY"] = 0.0;
      Process::environment.globals["SCSN-MP2 SAME-SPIN CORRELATION ENERGY"] = moinfo.emp2_ss*1.76;
      Process::environment.globals["SCSN-MP2 CORRELATION ENERGY"] = moinfo.emp2_ss*1.76;
      Process::environment.globals["SCSN-MP2 TOTAL ENERGY"] = moinfo.eref + moinfo.emp2_ss*1.76;
    }

    fprintf(outfile, "\n\tOpposite-spin MP2 correlation energy  = %20.15f\n", moinfo.emp2_os);
    fprintf(outfile, "\tSame-spin MP2 correlation energy      = %20.15f\n", moinfo.emp2_ss);
    fprintf(outfile, "\tMP2 correlation energy                = %20.15f\n", moinfo.emp2);
    fprintf(outfile, "      * MP2 total energy                      = %20.15f\n", moinfo.eref + moinfo.emp2);

    Process::environment.globals["MP2 OPPOSITE-SPIN CORRELATION ENERGY"] = moinfo.emp2_os;
    Process::environment.globals["MP2 SAME-SPIN CORRELATION ENERGY"] = moinfo.emp2_ss;
    Process::environment.globals["MP2 CORRELATION ENERGY"] = moinfo.emp2;
    Process::environment.globals["MP2 TOTAL ENERGY"] = moinfo.eref + moinfo.emp2;
  }
  if( params.wfn == "CC3"  || params.wfn == "EOM_CC3" ) {
    fprintf(outfile, "\tCC3 correlation energy     = %20.15f\n", moinfo.ecc);
    fprintf(outfile, "      * CC3 total energy           = %20.15f\n", moinfo.eref + moinfo.ecc);

    Process::environment.globals["CC3 CORRELATION ENERGY"] = moinfo.ecc;
    Process::environment.globals["CC3 TOTAL ENERGY"] = moinfo.eref + moinfo.ecc;
  }
  else if( params.wfn == "CC2" || params.wfn == "EOM_CC2" )  {
    fprintf(outfile, "\tCC2 correlation energy     = %20.15f\n", moinfo.ecc);
    fprintf(outfile, "      * CC2 total energy           = %20.15f\n", moinfo.eref + moinfo.ecc);

    Process::environment.globals["CC2 CORRELATION ENERGY"] = moinfo.ecc;
    Process::environment.globals["CC2 TOTAL ENERGY"] = moinfo.eref + moinfo.ecc;

    if(params.local && local.weakp == "MP2" )
      fprintf(outfile, "      * LCC2 (+LMP2) total energy  = %20.15f\n",
          moinfo.eref + moinfo.ecc + local.weak_pair_energy);
      Process::environment.globals["LCC2 (+LMP2) TOTAL ENERGY"] = 
          moinfo.eref + moinfo.ecc + local.weak_pair_energy;
  }
  else {
    if (params.scscc) {
      fprintf(outfile, "\n\tOS SCS-CCSD correlation energy      = %20.15f\n", moinfo.ecc_os*params.scscc_scale_os);
      fprintf(outfile, "\tSS SCS-CCSD correlation energy      = %20.15f\n", moinfo.ecc_ss*params.scscc_scale_ss);
      fprintf(outfile, "\tSCS-CCSD correlation energy            = %20.15f\n", moinfo.ecc_os*params.scscc_scale_os
            + moinfo.ecc_ss*params.scscc_scale_ss);
      fprintf(outfile, "      * SCS-CCSD total energy                  = %20.15f\n", moinfo.eref
            + moinfo.ecc_os*params.scscc_scale_os + moinfo.ecc_ss*params.scscc_scale_ss);

      // LAB TODO  reconsider variable names for ss/os cc
      Process::environment.globals["SCS-CCSD OPPOSITE-SPIN CORRELATION ENERGY"] = moinfo.ecc_os*params.scscc_scale_os;
      Process::environment.globals["SCS-CCSD SAME-SPIN CORRELATION ENERGY"] = moinfo.ecc_ss*params.scscc_scale_ss;
      Process::environment.globals["SCS-CCSD CORRELATION ENERGY"] = moinfo.ecc_os*params.scscc_scale_os + 
            moinfo.ecc_ss*params.scscc_scale_ss;
      Process::environment.globals["SCS-CCSD TOTAL ENERGY"] = moinfo.eref +
            moinfo.ecc_os*params.scscc_scale_os + moinfo.ecc_ss*params.scscc_scale_ss;
    }

    fprintf(outfile, "\n\tOpposite-spin CCSD correlation energy = %20.15f\n", moinfo.ecc_os);
    fprintf(outfile, "\tSame-spin CCSD correlation energy     = %20.15f\n", moinfo.ecc_ss);
    fprintf(outfile, "\tCCSD correlation energy               = %20.15f\n", moinfo.ecc);
    fprintf(outfile, "      * CCSD total energy                     = %20.15f\n", moinfo.eref + moinfo.ecc);

    Process::environment.globals["CCSD OPPOSITE-SPIN CORRELATION ENERGY"] = moinfo.ecc_os;
    Process::environment.globals["CCSD SAME-SPIN CORRELATION ENERGY"] = moinfo.ecc_ss;
    Process::environment.globals["CCSD CORRELATION ENERGY"] = moinfo.ecc;
    Process::environment.globals["CCSD TOTAL ENERGY"] = moinfo.ecc + moinfo.eref;

    if(params.local && local.weakp == "MP2" )
      fprintf(outfile, "      * LCCSD (+LMP2) total energy = %20.15f\n",
          moinfo.eref + moinfo.ecc + local.weak_pair_energy);
      Process::environment.globals["LCCSD (+LMP2) TOTAL ENERGY"] = 
          moinfo.eref + moinfo.ecc + local.weak_pair_energy;
  }
  fprintf(outfile, "\n");

  /* Write total energy to the checkpoint file */
  chkpt_init(PSIO_OPEN_OLD);
  chkpt_wt_etot(moinfo.ecc+moinfo.eref);
  chkpt_close();

  /* Write pertinent data to energy.dat for Dr. Yamaguchi */
//  if( params.wfn == "CCSD" || params.wfn == "BCCD" ) {
//
//    chkpt_init(PSIO_OPEN_OLD);
//    natom = chkpt_rd_natom();
//    geom = chkpt_rd_geom();
//    zvals = chkpt_rd_zvals();
//    chkpt_close();
//
//    ffile(&efile, "energy.dat",1);
//    fprintf(efile, "*\n");
//    for(i=0; i < natom; i++)
//      fprintf(efile, " %4d   %5.2f     %13.10f    %13.10f    %13.10f\n",
//          i+1, zvals[i], geom[i][0], geom[i][1], geom[i][2]);
//    free_block(geom);  free(zvals);
//    fprintf(efile, "SCF(30)   %22.12f\n", moinfo.escf);
//    fprintf(efile, "REF(100)  %22.12f\n", moinfo.eref);
//    if( params.wfn == "CCSD" )
//      fprintf(efile, "CCSD      %22.12f\n", (moinfo.ecc+moinfo.eref));
//    else if( params.wfn == "BCCD" )
//      fprintf(efile, "BCCD      %22.12f\n", (moinfo.ecc+moinfo.eref));
//    fclose(efile);
//  }

  /* Generate the spin-adapted RHF amplitudes for later codes */
  if(params.ref == 0) spinad_amps();

  /* Compute pair energies */
  if(params.print_pair_energies) {
    pair_energies(&ecc_aa, &ecc_ab);
    print_pair_energies(emp2_aa, emp2_ab, ecc_aa, ecc_ab);
  }

  if( (params.wfn == "CC3" || params.wfn == "EOM_CC3" )
      && (params.dertype == 1 || params.dertype == 3) && params.ref == 0) {
    params.ref = 1;
    /* generate the ROHF versions of the He^T1 intermediates */
    cc3_Wmnij();
    cc3_Wmbij();
    cc3_Wmnie();
    cc3_Wamef();
    cc3_Wabei();
    params.ref == 0;
  }

  if(params.local) {
    /*    local_print_T1_norm(); */
    local_done();
  }

  if(params.brueckner) brueckner_done = rotate();

  if( params.aobasis != "NONE" ) dpd_close(1);
  dpd_close(0);

  if(params.ref == 2) cachedone_uhf(cachelist);
  else cachedone_rhf(cachelist);
  free(cachefiles);

  cleanup();

#ifdef TIME_CCENERGY
  timer_off("CCEnergy");
#endif

  Process::environment.globals["CURRENT ENERGY"] = moinfo.ecc+moinfo.eref;

  exit_io();
  if(params.brueckner && brueckner_done)
      throw FeatureNotImplemented("CCENERGY", "Brueckner end loop", __FILE__, __LINE__);
  else
      return Success;
}

}} //namespace psi::ccenergy


namespace psi { namespace ccenergy {

void init_io()
{
//  int i, num_unparsed;
//  char *argv_unparsed[100];


  params.just_energy = 0;
  params.just_residuals = 0;
//  for (i=1, num_unparsed=0; i<argc; ++i) {
//    if (!strcmp(argv[i],"--just_energy")) {
//      /* just read T's on disk and compute energy */
//      params.just_energy = 1;
//    }
//    else if (!strcmp(argv[i],"--just_residuals")) {
//      params.just_residuals = 1;
//      /* just read T's on disk and compute residual matrix elements */
//    }
//    else {
//      argv_unparsed[num_unparsed++] = argv[i];
//    }
//  }

  tstart();

  for(int i=CC_MIN; i <= CC_MAX; i++) psio_open(i,1);
}

void title(void)
{
  fprintf(outfile, "\t\t\t**************************\n");
  fprintf(outfile, "\t\t\t*                        *\n");
  fprintf(outfile, "\t\t\t*        CCENERGY        *\n");
  fprintf(outfile, "\t\t\t*                        *\n");
  fprintf(outfile, "\t\t\t**************************\n");
}

void exit_io(void)
{
  int i;
  for(i=CC_MIN; i < CC_TMP; i++) psio_close(i,1);
  for(i=CC_TMP; i <= CC_TMP11; i++) psio_close(i,0); /* delete CC_TMP files */
  for(i=CC_TMP11+1; i <= CC_MAX; i++) psio_close(i,1);
  tstop();

}

void memchk(void)
{
  pid_t mypid;
  FILE *memdat;
  char comm[80];

  mypid = getpid();

  memdat = freopen("output.dat","a",stdout);

  sprintf(comm, "grep \"VmSize\" /proc/%d/status", (int) mypid);

  system(comm);
}

void init_ioff(void)
{
  int i;
  ioff = init_int_array(IOFF_MAX);
  ioff[0] = 0;
  for(i=1; i < IOFF_MAX; i++) ioff[i] = ioff[i-1] + i;
}


void checkpoint(void)
{
  int i;

  for(i=CC_MIN; i <= CC_MAX; i++) psio_close(i,1);
  for(i=CC_MIN; i <= CC_MAX; i++) psio_open(i,1);
}

/* just use T's on disk and don't iterate */
void one_step(void) {
  dpdfile2 t1;
    dpdbuf4 t2;
    double tval;

    moinfo.ecc = energy();
    fprintf(outfile,"\n\tValues computed from T amplitudes on disk.\n");
    fprintf(outfile,"Reference expectation value computed: %20.15lf\n", moinfo.ecc);
    psio_write_entry(CC_HBAR, "Reference expectation value", (char *) &(moinfo.ecc), sizeof(double));

  if (params.just_residuals) {
      Fme_build(); Fae_build(); Fmi_build();
        t1_build();
     Wmbej_build();
        Z_build();
        Wmnij_build();
        t2_build();
        if ( (params.ref == 0) || (params.ref == 1) ) {
          dpd_file2_init(&t1, CC_OEI, 0, 0, 1, "New tIA");
          dpd_file2_copy(&t1, CC_OEI, "FAI residual");
          dpd_file2_close(&t1);
          dpd_file2_init(&t1, CC_OEI, 0, 0, 1, "FAI residual");
            tval = dpd_file2_dot_self(&t1);
          dpd_file2_close(&t1);
            fprintf(outfile,"\tNorm squared of <Phi_I^A|Hbar|0> = %20.15lf\n",tval);
        }
      if (params.ref == 1) {
          dpd_file2_init(&t1, CC_OEI, 0, 0, 1, "New tia");
          dpd_file2_copy(&t1, CC_OEI, "Fai residual");
          dpd_file2_close(&t1);
        }
        else if (params.ref == 2) {
          dpd_file2_init(&t1, CC_OEI, 0, 2, 3, "New tia");
          dpd_file2_copy(&t1, CC_OEI, "Fai residual");
          dpd_file2_close(&t1);
        }
    if (params.ref == 0) {
          dpd_buf4_init(&t2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "New tIjAb");
          dpd_buf4_copy(&t2, CC_HBAR, "WAbIj residual");
          dpd_buf4_close(&t2);
            dpd_buf4_init(&t2, CC_HBAR, 0, 0, 5, 0, 5, 0, "WAbIj residual");
            tval = dpd_buf4_dot_self(&t2);
            fprintf(outfile,"\tNorm squared of <Phi^Ij_Ab|Hbar|0>: %20.15lf\n",tval);
            dpd_buf4_close(&t2);
    }
    else if (params.ref == 1) {
      dpd_buf4_init(&t2, CC_TAMPS, 0, 2, 7, 2, 7, 0, "New tIJAB");
      dpd_buf4_copy(&t2, CC_HBAR, "WABIJ residual");
      dpd_buf4_close(&t2);
      dpd_buf4_init(&t2, CC_TAMPS, 0, 2, 7, 2, 7, 0, "New tijab");
          dpd_buf4_copy(&t2, CC_HBAR, "Wabij residual");
      dpd_buf4_close(&t2);
      dpd_buf4_init(&t2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "New tIjAb");
          dpd_buf4_copy(&t2, CC_HBAR, "WAbIj residual");
      dpd_buf4_close(&t2);
    }
    else if(params.ref ==2) {
      dpd_buf4_init(&t2, CC_TAMPS, 0, 2, 7, 2, 7, 0, "New tIJAB");
      dpd_buf4_copy(&t2, CC_HBAR, "WABIJ residual");
      dpd_buf4_close(&t2);
      dpd_buf4_init(&t2, CC_TAMPS, 0, 12, 17, 12, 17, 0, "New tijab");
          dpd_buf4_copy(&t2, CC_HBAR, "Wabij residual");
      dpd_buf4_close(&t2);
      dpd_buf4_init(&t2, CC_TAMPS, 0, 22, 28, 22, 28, 0, "New tIjAb");
          dpd_buf4_copy(&t2, CC_HBAR, "WAbIj residual");
      dpd_buf4_close(&t2);
    }
    }
    return;
}

}} // namespace psi::ccenergy
