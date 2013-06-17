/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 *@END LICENSE
 */

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
int converged(double);
double diagnostic(void);
double d1diag(void);
double new_d1diag(void);
double d2diag(void);
void exit_io(void);
void cleanup(void);
void update(void);
void diis(int iter);
void ccdump(void);
int **cacheprep_uhf(int level, int *cachefiles);
int **cacheprep_rhf(int level, int *cachefiles);
void cachedone_rhf(int **cachelist);
void cachedone_uhf(int **cachelist);
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

    copy(reference_wavefunction_);

//    nso_        = reference_wavefunction_->nso();
//    nirrep_     = reference_wavefunction_->nirrep();
//    nmo_        = reference_wavefunction_->nmo();
//    for(int h = 0; h < nirrep_; ++h){
//        soccpi_[h] = reference_wavefunction_->soccpi()[h];
//        doccpi_[h] = reference_wavefunction_->doccpi()[h];
//        frzcpi_[h] = reference_wavefunction_->frzcpi()[h];
//        frzvpi_[h] = reference_wavefunction_->frzvpi()[h];
//        nmopi_[h]  = reference_wavefunction_->nmopi()[h];
//        nsopi_[h]  = reference_wavefunction_->nsopi()[h];
//    }
}

double CCEnergyWavefunction::compute_energy()
{
    energy_ = 0.0;
    PsiReturnType ccsd_return;
    if ((ccsd_return = psi::ccenergy::ccenergy(options_)) == Success) {
        // Get the total energy of the CCSD wavefunction
        energy_ = Process::environment.globals["CURRENT ENERGY"];
    }

    if ((options_.get_str("WFN") == "CCSD_T")) {
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

    std::vector<int*> spaces;
    spaces.push_back(moinfo.aoccpi);
    spaces.push_back(moinfo.aocc_sym);
    spaces.push_back(moinfo.avirtpi);
    spaces.push_back(moinfo.avir_sym);
    spaces.push_back(moinfo.boccpi);
    spaces.push_back(moinfo.bocc_sym);
    spaces.push_back(moinfo.bvirtpi);
    spaces.push_back(moinfo.bvir_sym);
    dpd_list[0] = boost::shared_ptr<DPD>( new DPD(0, moinfo.nirreps, params.memory, 0, cachefiles,
         cachelist, NULL, 4, spaces));
    dpd_set_default(0);

    if( params.aobasis != "NONE" ) { /* Set up new DPD's for AO-basis algorithm */
        std::vector<int*> aospaces;
        aospaces.push_back(moinfo.aoccpi);
        aospaces.push_back(moinfo.aocc_sym);
        aospaces.push_back(moinfo.sopi);
        aospaces.push_back(moinfo.sosym);
        aospaces.push_back(moinfo.boccpi);
        aospaces.push_back(moinfo.bocc_sym);
        aospaces.push_back(moinfo.sopi);
        aospaces.push_back(moinfo.sosym);
      dpd_list[1] = boost::shared_ptr<DPD>(new DPD(1, moinfo.nirreps, params.memory, 0, cachefiles,
                                                    cachelist, NULL, 4, aospaces));
      dpd_set_default(0);
    }

  }
  else { /** RHF or ROHF **/
    cachelist = cacheprep_rhf(params.cachelev, cachefiles);

    priority = priority_list();
    std::vector<int*> spaces;
    spaces.push_back(moinfo.occpi);
    spaces.push_back(moinfo.occ_sym);
    spaces.push_back(moinfo.virtpi);
    spaces.push_back(moinfo.vir_sym);

    dpd_list[0] = boost::shared_ptr<DPD>(new DPD(0, moinfo.nirreps, params.memory, params.cachetype, cachefiles,
         cachelist, priority, 2, spaces));
    dpd_set_default(0);


    if( params.aobasis != "NONE") { /* Set up new DPD for AO-basis algorithm */
            std::vector<int*> aospaces;
            aospaces.push_back(moinfo.occpi);
            aospaces.push_back(moinfo.occ_sym);
            aospaces.push_back(moinfo.sopi);
            aospaces.push_back(moinfo.sosym);
      dpd_list[1] = boost::shared_ptr<DPD>(new DPD(1, moinfo.nirreps, params.memory, 0, cachefiles, cachelist, NULL,
               2, aospaces));
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
    psio_write_entry(PSIF_CC_INFO, "MP2 Energy", (char *) &(moinfo.emp2),sizeof(double));
    Process::environment.globals["MP2 CORRELATION ENERGY"] = moinfo.emp2;
    Process::environment.globals["MP2 TOTAL ENERGY"] = moinfo.emp2 + moinfo.eref;
  }

  if(params.print_mp2_amps) amp_write();

  tau_build();
  taut_build();
  fprintf(outfile, "\t            Solving CC Amplitude Equations\n");
  fprintf(outfile, "\t            ------------------------------\n");
  fprintf(outfile, "  Iter             Energy              RMS        T1Diag      D1Diag    New D1Diag    D2Diag\n");
  fprintf(outfile, "  ----     ---------------------    ---------   ----------  ----------  ----------   --------\n");
  moinfo.ecc = energy();
  pair_energies(&emp2_aa, &emp2_ab);
  double last_energy;

  moinfo.t1diag = diagnostic();
  moinfo.d1diag = d1diag();
  moinfo.new_d1diag = new_d1diag();
  fflush(outfile);
  moinfo.d2diag = d2diag();
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

    if(converged(last_energy - moinfo.ecc)) {
      done = 1;

      tsave();
      tau_build(); taut_build();
      last_energy = moinfo.ecc;
      moinfo.ecc = energy();
      moinfo.t1diag = diagnostic();
      moinfo.d1diag = d1diag();
      moinfo.new_d1diag = new_d1diag();
      moinfo.d2diag = d2diag();
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
    last_energy = moinfo.ecc;
    moinfo.ecc = energy();
    moinfo.t1diag = diagnostic();
    moinfo.d1diag = d1diag();
    moinfo.new_d1diag = new_d1diag();
    moinfo.d2diag = d2diag();
    update();
    checkpoint();
  }  // end loop over iterations
  fprintf(outfile, "\n");
  if(!done) {
    fprintf(outfile, "\t ** Wave function not converged to %2.1e ** \n",
        params.convergence);
    fflush(outfile);
//    if( params.aobasis != "NONE" ) dpd_->close(1);
//    dpd_->close(0);
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
//    params.ref == 0;
  }

  if(params.local) {
    /*    local_print_T1_norm(); */
    local_done();
  }

  if(params.brueckner)
    Process::environment.globals["BRUECKNER CONVERGED"] = rotate();

//  if( params.aobasis != "NONE" ) dpd_->close(1);
//  dpd_->close(0);

  if(params.ref == 2) cachedone_uhf(cachelist);
  else cachedone_rhf(cachelist);
  free(cachefiles);

  cleanup();

#ifdef TIME_CCENERGY
  timer_off("CCEnergy");
#endif

  Process::environment.globals["CURRENT ENERGY"] = moinfo.ecc+moinfo.eref;
  Process::environment.globals["CURRENT CORRELATION ENERGY"] = moinfo.ecc;
  //Process::environment.globals["CC TOTAL ENERGY"] = moinfo.ecc+moinfo.eref;
  //Process::environment.globals["CC CORRELATION ENERGY"] = moinfo.ecc;

  exit_io();
//  if(params.brueckner && brueckner_done)
 //     throw FeatureNotImplemented("CCENERGY", "Brueckner end loop", __FILE__, __LINE__);
  //else
  return Success;
}

}} //namespace psi::ccenergy


namespace psi { namespace ccenergy {

void init_io()
{
  params.just_energy = 0;
  params.just_residuals = 0;
  tstart();
  for(int i =PSIF_CC_MIN; i <= PSIF_CC_MAX; i++) psio_open(i,1);
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
  for(i=PSIF_CC_MIN; i < PSIF_CC_TMP; i++) psio_close(i,1);
  for(i=PSIF_CC_TMP; i <= PSIF_CC_TMP11; i++) psio_close(i,0); /* delete CC_TMP files */
  for(i=PSIF_CC_TMP11+1; i <= PSIF_CC_MAX; i++) psio_close(i,1);
  tstop();

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

  for(i=PSIF_CC_MIN; i <= PSIF_CC_MAX; i++) psio_close(i,1);
  for(i=PSIF_CC_MIN; i <= PSIF_CC_MAX; i++) psio_open(i,1);
}

/* just use T's on disk and don't iterate */
void one_step(void) {
  dpdfile2 t1;
    dpdbuf4 t2;
    double tval;

    moinfo.ecc = energy();
    fprintf(outfile,"\n\tValues computed from T amplitudes on disk.\n");
    fprintf(outfile,"Reference expectation value computed: %20.15lf\n", moinfo.ecc);
    psio_write_entry(PSIF_CC_HBAR, "Reference expectation value", (char *) &(moinfo.ecc), sizeof(double));

  if (params.just_residuals) {
      Fme_build(); Fae_build(); Fmi_build();
        t1_build();
     Wmbej_build();
        Z_build();
        Wmnij_build();
        t2_build();
        if ( (params.ref == 0) || (params.ref == 1) ) {
          dpd_->file2_init(&t1, PSIF_CC_OEI, 0, 0, 1, "New tIA");
          dpd_->file2_copy(&t1, PSIF_CC_OEI, "FAI residual");
          dpd_->file2_close(&t1);
          dpd_->file2_init(&t1, PSIF_CC_OEI, 0, 0, 1, "FAI residual");
            tval = dpd_->file2_dot_self(&t1);
          dpd_->file2_close(&t1);
            fprintf(outfile,"\tNorm squared of <Phi_I^A|Hbar|0> = %20.15lf\n",tval);
        }
      if (params.ref == 1) {
          dpd_->file2_init(&t1, PSIF_CC_OEI, 0, 0, 1, "New tia");
          dpd_->file2_copy(&t1, PSIF_CC_OEI, "Fai residual");
          dpd_->file2_close(&t1);
        }
        else if (params.ref == 2) {
          dpd_->file2_init(&t1, PSIF_CC_OEI, 0, 2, 3, "New tia");
          dpd_->file2_copy(&t1, PSIF_CC_OEI, "Fai residual");
          dpd_->file2_close(&t1);
        }
    if (params.ref == 0) {
          dpd_->buf4_init(&t2, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "New tIjAb");
          dpd_->buf4_copy(&t2, PSIF_CC_HBAR, "WAbIj residual");
          dpd_->buf4_close(&t2);
            dpd_->buf4_init(&t2, PSIF_CC_HBAR, 0, 0, 5, 0, 5, 0, "WAbIj residual");
            tval = dpd_->buf4_dot_self(&t2);
            fprintf(outfile,"\tNorm squared of <Phi^Ij_Ab|Hbar|0>: %20.15lf\n",tval);
            dpd_->buf4_close(&t2);
    }
    else if (params.ref == 1) {
      dpd_->buf4_init(&t2, PSIF_CC_TAMPS, 0, 2, 7, 2, 7, 0, "New tIJAB");
      dpd_->buf4_copy(&t2, PSIF_CC_HBAR, "WABIJ residual");
      dpd_->buf4_close(&t2);
      dpd_->buf4_init(&t2, PSIF_CC_TAMPS, 0, 2, 7, 2, 7, 0, "New tijab");
          dpd_->buf4_copy(&t2, PSIF_CC_HBAR, "Wabij residual");
      dpd_->buf4_close(&t2);
      dpd_->buf4_init(&t2, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "New tIjAb");
          dpd_->buf4_copy(&t2, PSIF_CC_HBAR, "WAbIj residual");
      dpd_->buf4_close(&t2);
    }
    else if(params.ref ==2) {
      dpd_->buf4_init(&t2, PSIF_CC_TAMPS, 0, 2, 7, 2, 7, 0, "New tIJAB");
      dpd_->buf4_copy(&t2, PSIF_CC_HBAR, "WABIJ residual");
      dpd_->buf4_close(&t2);
      dpd_->buf4_init(&t2, PSIF_CC_TAMPS, 0, 12, 17, 12, 17, 0, "New tijab");
          dpd_->buf4_copy(&t2, PSIF_CC_HBAR, "Wabij residual");
      dpd_->buf4_close(&t2);
      dpd_->buf4_init(&t2, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "New tIjAb");
          dpd_->buf4_copy(&t2, PSIF_CC_HBAR, "WAbIj residual");
      dpd_->buf4_close(&t2);
    }
    }
    return;
}

}} // namespace psi::ccenergy
