/*! \defgroup MP2 mp2: Canonical Evaluation of MP2 Energy and Gradients */

/*!
** \file
** \ingroup MP2
** \brief Canonical evaluation of MP2 energy and gradients
*/

#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <libciomr/libciomr.h>
#include <libdpd/dpd.h>
#include <libchkpt/chkpt.h>
#include <libqt/qt.h>
#include <libiwl/iwl.h>
#include <physconst.h>
#include <psifiles.h>
#include <libmints/wavefunction.h>
#include "mp2wave.h"
#include "globals.h"

namespace psi{ namespace mp2{

void init_io();
void title(void);
void get_moinfo(void);
void get_params(Options &options);
void init_ioff(void);
int **cacheprep_rhf(int level, int *cachefiles);
int **cacheprep_uhf(int level, int *cachefiles);
void cachedone_rhf(int **cachelist);
double energy(void);
void amps(void);
void sort_amps(void);
void opdm(void);
void twopdm(void);
void lag(void);
void build_X(void);
void build_A(void);
void Zvector(void);
void relax_I(void);
void relax_opdm(void);
void sort_opdm(void);
void sort_I(void);
void fold(void);
void deanti(void);
void write_data(void);
void check_energy(int);
void sort_twopdm(void);
void cleanup(void);
void exit_io(void);

PsiReturnType
mp2(Options & options)
{
  using namespace psi::mp2;

  int *cachefiles;
  int **cachelist;

  init_io();
  title();

  get_moinfo();
  get_params(options);
  init_ioff();

  cachefiles = init_int_array(PSIO_MAXUNIT);

  if(params.ref == 2) { /** UHF **/
    cachelist = cacheprep_uhf(params.cachelev,cachefiles);
    dpd_init(0,mo.nirreps,params.memory,params.cachetype,cachefiles,cachelist,
             NULL,4,mo.aoccpi,mo.aocc_sym,mo.avirpi,mo.avir_sym,mo.boccpi,
             mo.bocc_sym,mo.bvirpi,mo.bvir_sym);
  }
  else { /** RHF or ROHF **/
    cachelist = cacheprep_rhf(params.cachelev,cachefiles);
    dpd_init(0,mo.nirreps,params.memory,params.cachetype,cachefiles,cachelist,
         NULL,2,mo.occpi,mo.occ_sym,mo.virpi,mo.vir_sym);
  }

  amps();

  mo.Emp2 = energy();

  fprintf(outfile,"\n");
  fprintf(outfile,"\tScaled_OS correlation energy      = %20.15f\n",mo.escsmp2_os);
  fprintf(outfile,"\tScaled_SS correlation energy      = %20.15f\n",mo.escsmp2_ss);
  fprintf(outfile,"\tSCS-MP2 correlation energy        = %20.15f\n",mo.escsmp2_os+mo.escsmp2_ss);
  fprintf(outfile,"      * SCS-MP2 total energy              = %20.15f\n",mo.Escf+mo.escsmp2_os+mo.escsmp2_ss);

  fprintf(outfile,"\n\tOpposite-spin correlation energy  = %20.15f\n",mo.emp2_os);
  fprintf(outfile,"\tSame-spin correlation energy      = %20.15f\n",mo.emp2_ss);
  fprintf(outfile,"\tMP2 correlation energy            = %20.15f\n",mo.Emp2);
  fprintf(outfile,"      * MP2 total energy                  = %20.15f\n\n",mo.Escf + mo.Emp2);
  fflush(outfile);

  Process::environment.globals["CURRENT ENERGY"] = mo.Escf + mo.Emp2;
  Process::environment.globals["SCS-MP2 CORRELATION ENERGY"] = mo.escsmp2_os+mo.escsmp2_ss;
  Process::environment.globals["SAME SPIN MP2 CORRELATION ENERGY"] = mo.escsmp2_ss;
  Process::environment.globals["OPPOSITE SPIN MP2 CORRELATION ENERGY"] = mo.escsmp2_os;
  Process::environment.globals["MP2 CORRELATION ENERGY"] = mo.Escf + mo.Emp2;

  chkpt_init(PSIO_OPEN_OLD);
  // Save MP2 contribution to Chkpt
  chkpt_wt_emp2(mo.Emp2);
  chkpt_wt_etot(mo.Escf+mo.Emp2);
  chkpt_close();

  if(params.opdm) {
    sort_amps();
    opdm();
    if(params.relax_opdm) {
      lag();
      build_A();
      Zvector();
    }
    sort_opdm();
    //dipole();
  }

  if(params.gradient) {
    sort_amps();
    opdm();
    twopdm();
    lag();
    build_X();
    build_A();
    Zvector();
    relax_I();
    relax_opdm();
    sort_I();
    sort_opdm();
    fold();
    deanti();
    write_data();
  }

  dpd_close(0);

  cleanup();

  exit_io();

  return Success;
}


void init_io()
{
//  int i=0;
//  int num_extra_args=0;
//  char **extra_args;

//  extra_args = (char **) malloc(argc*sizeof(char *));
  params.opdm = 0;

//  for(i=1; i<argc; i++) {
//    if(strcmp(argv[i], "--opdm") == 0) {
//      params.opdm = 1;
//    }
//    else {
//      extra_args[num_extra_args++] = argv[i];
//    }
//  }

  tstart();

  for(int i=CC_MIN; i <= CC_MAX; i++)
    psio_open(i,1);

//  free(extra_args);
}

void title(void)
{
  fprintf(outfile, "\t\t\t*************************\n");
  fprintf(outfile, "\t\t\t*                       *\n");
  fprintf(outfile, "\t\t\t*          MP2          *\n");
  fprintf(outfile, "\t\t\t*                       *\n");
  fprintf(outfile, "\t\t\t*************************\n");
  fflush(outfile);
}

void init_ioff(void)
{
  int i;

  ioff = init_int_array(MAXIOFF);
  ioff[0] = 0;
  for(i=1; i < MAXIOFF; i++) {
    ioff[i] = ioff[i-1] + i;
  }

}

void cleanup(void)
{
  int i;

  free(mo.doccpi);
  free(mo.soccpi);
  free(mo.mopi);
  free(mo.fzdoccpi);
  free(mo.fzvirtpi);
  for(i=0; i < mo.nirreps; i++)
    free(mo.irreplabels[i]);
  free(mo.irreplabels);
  free(ioff);

  if(params.ref == 2) {
    free(mo.aoccpi);
    free(mo.boccpi);
    free(mo.avirpi);
    free(mo.bvirpi);
    free(mo.aocc_sym);
    free(mo.bocc_sym);
    free(mo.avir_sym);
    free(mo.bvir_sym);
    free(mo.aocc_off);
    free(mo.bocc_off);
    free(mo.avir_off);
    free(mo.bvir_off);
    free(mo.qt_aocc);
    free(mo.qt_bocc);
    free(mo.qt_avir);
    free(mo.qt_bvir);
  }
  else {
    free(mo.occpi);
    free(mo.virpi);
    free(mo.occ_sym);
    free(mo.vir_sym);
    free(mo.occ_off);
    free(mo.vir_off);
    free(mo.qt_occ);
    free(mo.qt_vir);
  }
}

void exit_io(void)
{
  int i;
  for(int i=CC_MIN; i <= CC_MAX; i++)
    psio_close(i,0);

  tstop();
}

MP2Wavefunction::MP2Wavefunction(boost::shared_ptr<Wavefunction> reference_wavefunction, Options &options)
    : Wavefunction(options, _default_psio_lib_)
{
    set_reference_wavefunction(reference_wavefunction);
    init();
}

MP2Wavefunction::~MP2Wavefunction()
{

}

void MP2Wavefunction::init()
{
    // Wavefunction creates a chkpt object for your, but we're not going to use it.
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

double MP2Wavefunction::compute_energy()
{
    energy_ = 0.0;
    PsiReturnType mp2_return;
    if ((mp2_return = psi::mp2::mp2(options_)) == Success) {
        energy_ = Process::environment.globals["CURRENT ENERGY"];
    }

    return energy_;
}
}} /* End namespaces */
