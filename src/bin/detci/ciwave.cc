#include <libmints/wavefunction.h>
#include <libmints/matrix.h>
#include <libiwl/iwl.h>
#include <psifiles.h>
#include "ciwave.h"
#include "structs.h"
#define EXTERN
#include "globals.h"
#include "globaldefs.h"


namespace psi { namespace detci {

CIWavefunction::CIWavefunction(boost::shared_ptr<Wavefunction> reference_wavefunction,
                               Options &options)
    : Wavefunction(options, reference_wavefunction->psio())
{
    set_reference_wavefunction(reference_wavefunction);

    // TODO-CDS:
    // The CC codes destroy the checkpoint object created by Wavefunction.
    // We'd like to be able to do the same here.  Need to convert everything
    // such that we don't explicitly need checkpoint
    // Destroy it. Otherwise we will see a "file already open" error.
    chkpt_.reset();
    set_opdm(false, true);
    set_orbitals(false, true);

}

CIWavefunction::~CIWavefunction()
{

}

void CIWavefunction::common_init()
{

    // We're copying this stuff over like it's done in ccenergy, but
    // these Wavefunction member data are not actually used by the code
    // yet (Sept 2011 CDS).  Copying these only in case they're needed
    // by someone else who uses the CIWavefunction and expects them to
    // be available.

    // CDS-TODO: Note, some of these data should be updated to reflect what
    // the CI wavefunction is doing, not what the reference wavefunction
    // is doing, in case these are actually used elsewhere.
    nso_        = reference_wavefunction_->nso();
    nirrep_     = reference_wavefunction_->nirrep();
    nmo_        = reference_wavefunction_->nmo();

    nalpha_     = CalcInfo.num_alp;
    nbeta_      = CalcInfo.num_bet;

    // Per irrep data
    for(int h = 0; h < nirrep_; ++h){
        soccpi_[h] = CalcInfo.socc[h];
        doccpi_[h] = CalcInfo.docc[h];
        frzcpi_[h] = reference_wavefunction_->frzcpi()[h];
        frzvpi_[h] = reference_wavefunction_->frzvpi()[h];
        nmopi_[h]  = CalcInfo.orbs_per_irr[h];
        nsopi_[h]  = CalcInfo.so_per_irr[h];
    }
}

double CIWavefunction::compute_energy()
{
    energy_ = 0.0;
    // PsiReturnType ci_return;
    // if ((ci_return = detci(options_)) == Success) {
    //     // Get the total energy
    //     energy_ = Process::environment.globals["CURRENT ENERGY"];
    // }

    return energy_;
}

void CIWavefunction::set_opdm(bool print, bool erase)
{
    int npop = CalcInfo.num_ci_orbs + CalcInfo.num_fzc_orbs;
    SharedMatrix opdm_a(new Matrix("One-Particle Alpha Density Matrix", npop, npop));
    double *opdm_ap = opdm_a->pointer()[0];
  
    SharedMatrix opdm_b(new Matrix("One-Particle Beta Density Matrix", npop, npop));
    double *opdm_bp = opdm_b->pointer()[0];
  
    char opdm_key[80];
    int root = 1;
  
    psio_open(Parameters.opdm_file, PSIO_OPEN_OLD);
    /* if the user hasn't specified a root, just get "the" onepdm */
    if (!options_["FOLLOW_ROOT"].has_changed()) {
        psio_read_entry(Parameters.opdm_file, "MO-basis Alpha OPDM", (char *) opdm_ap,
                        npop*npop*sizeof(double));
        psio_read_entry(Parameters.opdm_file, "MO-basis Beta OPDM", (char *) opdm_bp,
                        npop*npop*sizeof(double));
    }
    else {
        root = options_.get_int("FOLLOW_ROOT");
        sprintf(opdm_key, "MO-basis Alpha OPDM Root %d", root);
        psio_read_entry(Parameters.opdm_file, opdm_key, (char *) opdm_ap,
                        npop*npop*sizeof(double));
        sprintf(opdm_key, "MO-basis Beta OPDM Root %d", root);
        psio_read_entry(Parameters.opdm_file, opdm_key, (char *) opdm_bp,
                        npop*npop*sizeof(double));
    }
  
    psio_close(Parameters.opdm_file, erase ? 0 : 1);
  
    if (print){
        opdm_a->print();
        opdm_b->print();
    }

    Da_ = opdm_a;
    Db_ = opdm_b;
}

void CIWavefunction::set_orbitals(bool print, bool erase)
{
    int h, ir_orbs;
    char orb_key[80];

    SharedMatrix orbitals(new Matrix("Orbitals", CalcInfo.nirreps, 
                          CalcInfo.orbs_per_irr, CalcInfo.orbs_per_irr));
    double *orbp;
    
  
    psio_open(PSIF_DETCAS, PSIO_OPEN_OLD);
    for (h=0; h<CalcInfo.nirreps; h++) {
      ir_orbs = CalcInfo.orbs_per_irr[h];
      sprintf(orb_key, "Orbs Irrep %2d", h);
      psio_read_entry(PSIF_DETCAS, orb_key, (char *) orbitals->pointer(h)[0],
                      ir_orbs*ir_orbs*sizeof(double));
    }
    psio_close(PSIF_DETCAS, erase ? 0 : 1);

    if (print) orbitals->print();
    Ca_ = orbitals;
     
}
// void CIWavefunction::finalize()
// {
// 
// }


}}
