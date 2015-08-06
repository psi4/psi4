#include <libmints/mints.h>
#include <psi4-dec.h>
#include <libiwl/iwl.h>
#include <libciomr/libciomr.h>
#include <psifiles.h>
#include <libqt/qt.h>

#include "globaldefs.h"
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

    common_init();
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

    get_mo_info();        /* read DOCC, SOCC, frozen, nmo, etc        */

    // Wavefunction frozen nomenclature is equivalent to dropped in detci.
    // In detci frozen means doubly occupied, but no orbital rotations.
    nirrep_     = reference_wavefunction_->nirrep();
    nso_        = reference_wavefunction_->nso();
    nmo_        = reference_wavefunction_->nmo();
    nalpha_     = CalcInfo.num_alp; // Total number of alpha electrons, including core
    nbeta_      = CalcInfo.num_bet;
    nfrzc_      = CalcInfo.num_drc_orbs;

    // Per irrep data
    for(int h = 0; h < nirrep_; ++h){
        frzcpi_[h] = CalcInfo.dropped_docc[h];
        doccpi_[h] = CalcInfo.docc[h];
        soccpi_[h] = CalcInfo.socc[h];
        frzvpi_[h] = CalcInfo.dropped_uocc[h];

        nmopi_[h]  = CalcInfo.orbs_per_irr[h];
        nsopi_[h]  = CalcInfo.so_per_irr[h];

    H_ = reference_wavefunction_->H();
    Ca_ = reference_wavefunction_->Ca()->clone();
    psio_ = reference_wavefunction_->psio();
    AO2SO_ = reference_wavefunction_->aotoso();
    molecule_ = reference_wavefunction_->molecule();

    name_ = "CIWavefunction";
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

void CIWavefunction::orbital_locations(const std::string& orbitals, int* start, int* end){
    if (orbitals == "FZC"){
      for (int h=0; h<nirrep_; h++){
        start[h] = 0;
        end[h] = CalcInfo.frozen_docc[h];
      }
    }
    else if (orbitals == "DOCC"){
      for (int h=0; h<nirrep_; h++){
        start[h] = CalcInfo.frozen_docc[h];
        end[h] = CalcInfo.dropped_docc[h];
      }
    }
    else if (orbitals == "DRC"){
      for (int h=0; h<nirrep_; h++){
        start[h] = 0;
        end[h] = CalcInfo.dropped_docc[h];
      }
    }
    else if (orbitals == "ACT"){
      for (int h=0; h<nirrep_; h++){
        start[h] = CalcInfo.dropped_docc[h];
        end[h] = nmopi_[h] - CalcInfo.dropped_uocc[h];
      }
    }
    else if (orbitals == "DRV"){
      for (int h=0; h<nirrep_; h++){
        start[h] = nmopi_[h] - CalcInfo.dropped_uocc[h];
        end[h] = nmopi_[h];
      }
    }
    else if (orbitals == "VIR"){
      for (int h=0; h<nirrep_; h++){
        start[h] = nmopi_[h] - CalcInfo.dropped_uocc[h];
        end[h] = nmopi_[h] - CalcInfo.frozen_uocc[h];
      }
    }
    else if (orbitals == "FZV"){
      for (int h=0; h<nirrep_; h++){
        start[h] = nmopi_[h] - CalcInfo.frozen_uocc[h];
        end[h] = nmopi_[h];
      }
    }
    else if (orbitals == "ALL"){
      for (int h=0; h<nirrep_; h++){
        start[h] = 0;
        end[h] = nmopi_[h];
      }
    }
    else if (orbitals == "ROT"){
      for (int h=0; h<nirrep_; h++){
        start[h] = CalcInfo.frozen_docc[h];
        end[h] = nmopi_[h] - CalcInfo.frozen_uocc[h];
      }
    }
    else{
        throw PSIEXCEPTION("DETCI: Orbital subset is not defined, should be FZC, DOCC, ACT, VIR, FZV, or ALL");
    }
}

SharedMatrix CIWavefunction::get_orbitals(const std::string& orbital_name)
{
    /// Figure out orbital positions
    int* start = new int[nirrep_];
    int* end = new int[nirrep_];

    orbital_locations(orbital_name, start, end);

    int* spread = new int[nirrep_];
    for (int h=0; h<nirrep_; h++){
      spread[h] = end[h] - start[h];
    }

    /// Fill desired orbitals
    SharedMatrix retC(new Matrix("C " + orbital_name, nirrep_, nsopi_, spread));
    for (int h = 0; h < nirrep_; h++) {
        for (int i = start[h], pos=0; i < end[h]; i++, pos++) {
            C_DCOPY(nsopi_[h], &Ca_->pointer(h)[0][i], nmopi_[h], &retC->pointer(h)[0][pos], spread[h]);
        }
    }

    /// Cleanup
    delete[] start;
    delete[] end;
    delete[] spread;

    return retC;
}

void CIWavefunction::set_orbitals(const std::string& orbital_name, SharedMatrix orbitals)
{
    /// Figure out orbital positions
    int* start = new int[nirrep_];
    int* end = new int[nirrep_];

    orbital_locations(orbital_name, start, end);

    int* spread = new int[nirrep_];
    for (int h=0; h<nirrep_; h++){
      spread[h] = end[h] - start[h];
    }

    /// Fill desired orbitals
    for (int h = 0; h < nirrep_; h++) {
        for (int i = start[h], pos=0; i < end[h]; i++, pos++) {
            C_DCOPY(nsopi_[h], &orbitals->pointer(h)[0][pos], spread[h], &Ca_->pointer(h)[0][i], nmopi_[h]);
        }
    }

    /// Cleanup
    delete[] start;
    delete[] end;
    delete[] spread;
}
void CIWavefunction::set_opdm(bool use_old_d)
{

  // No new opdm
  if(use_old_d){
    Da_ = reference_wavefunction_->Da();
    Db_ = reference_wavefunction_->Db();
    return;
  }

  int npop = CalcInfo.num_ci_orbs + CalcInfo.num_drc_orbs;
  SharedMatrix opdm_a(new Matrix("MO-basis Alpha OPDM", npop, npop));
  double *opdm_ap = opdm_a->pointer()[0];

  SharedMatrix opdm_b(new Matrix("MO-basis Beta OPDM", npop, npop));
  double *opdm_bp = opdm_b->pointer()[0];

  char opdm_key[80];
  psio_open(Parameters.opdm_file, PSIO_OPEN_OLD);
  if (!options_["FOLLOW_ROOT"].has_changed()) {
    psio_read_entry(Parameters.opdm_file, "MO-basis Alpha OPDM", (char *) opdm_ap,
                    npop*npop*sizeof(double));
    psio_read_entry(Parameters.opdm_file, "MO-basis Beta OPDM", (char *) opdm_bp,
                    npop*npop*sizeof(double));
  }
  else {
    int root = options_.get_int("FOLLOW_ROOT");
    sprintf(opdm_key, "MO-basis OPDM Alpha Root %d", root);
    psio_read_entry(Parameters.opdm_file, opdm_key, (char *) opdm_ap,
                    npop*npop*sizeof(double));
    sprintf(opdm_key, "MO-basis OPDM Beta Root %d", root);
    psio_read_entry(Parameters.opdm_file, opdm_key, (char *) opdm_bp,
                    npop*npop*sizeof(double));
  }
  psio_close(Parameters.opdm_file, 1);

  Da_ = opdm_a;
  Db_ = opdm_b;

}
SharedMatrix CIWavefunction::get_opdm(int root, int spin, bool transden)
{
    int npop = CalcInfo.num_ci_orbs + CalcInfo.num_drc_orbs;
    SharedMatrix opdm_a(new Matrix("One-Particle Alpha Density Matrix", npop, npop));
    double *opdm_ap = opdm_a->pointer()[0];

    SharedMatrix opdm_b(new Matrix("One-Particle Beta Density Matrix", npop, npop));
    double *opdm_bp = opdm_b->pointer()[0];

    char opdm_key[80];
    psio_open(Parameters.opdm_file, PSIO_OPEN_OLD);
    sprintf(opdm_key, "MO-basis Alpha OPDM %s Root %d", transden ? "TDM" : "OPDM", root);
    psio_read_entry(Parameters.opdm_file, opdm_key, (char *) opdm_ap,
                    npop*npop*sizeof(double));
    sprintf(opdm_key, "MO-basis Beta OPDM %s Root %d", transden ? "TDM" : "OPDM", root);
    psio_read_entry(Parameters.opdm_file, opdm_key, (char *) opdm_bp,
                        npop*npop*sizeof(double));
    psio_close(Parameters.opdm_file, 1);

    if (spin == 2){
        opdm_a->add(opdm_b);
        return opdm_a;
    }
    else if (spin == 1){
        return opdm_b;
    }
    else if (spin == 0){
        return opdm_a;
    }
    else {
       throw PSIEXCEPTION("DETCI: option spin is not recognizable value.");
    }
}
SharedMatrix CIWavefunction::get_active_opdm()
{
  SharedMatrix actOPDM(new Matrix("OPDM", nirrep_, CalcInfo.ci_orbs, CalcInfo.ci_orbs));
  double** OPDMa = Da_->pointer();
  double** OPDMb = Db_->pointer();

  int offset = 0;
  for (int h=0; h<nirrep_; h++){
    offset += CalcInfo.dropped_docc[h];
    if (!CalcInfo.ci_orbs[h]) continue;

    double* actp = actOPDM->pointer(h)[0];

    for (int i=0, target=0; i<CalcInfo.ci_orbs[h]; i++){
      int offi = offset + i;
      for (int j=0; j<CalcInfo.ci_orbs[h]; j++){
        actp[target++] = OPDMa[offi][offset + j] + OPDMb[offi][offset + j];
      }
    }
    offset += CalcInfo.ci_orbs[h];
  }
  return actOPDM;
}
void CIWavefunction::set_lag()
{
    SharedMatrix lag(new Matrix("MCSCF Lagrangian", nmo_, nmo_));
    double **lagp = lag->pointer();
    for (size_t i=0; i<nmo_; i++){
        for (size_t j=0; j<nmo_; j++){
            lagp[i][j] = MCSCF_CalcInfo.lag[i][j];
        }
    }
    Lagrangian_ = lag;
}
SharedVector CIWavefunction::get_tpdm(bool symmetrize, bool act_only)
{

  int sqnbf, ntri;
  int *ioff_lt, i;                    /* offsets for left (or right) indices */
  int p,q,r,s,smax,pq,qp,rs,sr,pqrs,qprs,pqsr,qpsr,target;
  struct iwlbuf TBuff;

  iwl_buf_init(&TBuff, Parameters.tpdm_file, 0.0, 1, 1);
  int npop = CalcInfo.num_ci_orbs + CalcInfo.num_drc_orbs;

  sqnbf = npop * npop;
  SharedVector tpdm(new Vector("TPDM", (sqnbf*(sqnbf+1))/2));
  double* tpdmp = tpdm->pointer();

  /* Construct the ioff_lt array (same here as ioff_rt) : different than
   * regular ioff because there is no perm symmetry between left indices
   * or right indices.
   */
  ioff_lt = init_int_array(nmo_);
  for (i=0; i<npop; i++) {
    ioff_lt[i] = i * npop;
  }

 iwl_buf_rd_all(&TBuff, tpdmp, ioff_lt, ioff_lt, 1, ioff,
                false, "outfile");

  iwl_buf_close(&TBuff, 1);
  free(ioff_lt);

  if(!symmetrize){
    return tpdm;
  }

  ntri = (npop * (npop + 1))/2;
  SharedVector symm_tpdm(new Vector("Symm TPDM", (ntri * (ntri + 1))/2));
  double* symm_tpdmp = symm_tpdm->pointer();

  for (p=0,target=0; p<npop; p++) {
    for (q=0; q<=p; q++) {
      for (r=0; r<=p; r++) {
        smax = (r==p) ? q+1 : r+1;
        for (s=0; s<smax; s++,target++) {

          pq = p * npop + q;
          qp = q * npop + p;
          rs = r * npop + s;
          sr = s * npop + r;
          pqrs = INDEX(pq,rs);
          qprs = INDEX(qp,rs);
          pqsr = INDEX(pq,sr);
          qpsr = INDEX(qp,sr);
          /* would be 0.25 but the formulae I used for the diag hessian
           * seem to define the TPDM with the 1/2 back outside */
          symm_tpdmp[target] = 0.5 * (tpdmp[pqrs] + tpdmp[qprs] +
                         tpdmp[pqsr] + tpdmp[qpsr]);

        }
      }
    }
  }
  tpdm.reset();
  return symm_tpdm;
}
SharedMatrix CIWavefunction::get_active_tpdm()
{
  int nci   = CalcInfo.num_ci_orbs;
  int ndocc = CalcInfo.num_drc_orbs;
  int npop  = CalcInfo.num_ci_orbs + CalcInfo.num_drc_orbs;
  SharedMatrix actTPDM(new Matrix("TPDM", nci*nci, nci*nci));
  double* actTPDMp = actTPDM->pointer()[0];

  SharedVector fullTPDM = get_tpdm();
  double* fullTPDMp = fullTPDM->pointer();

  int nci2 = nci * nci;
  int nci3 = nci * nci * nci;
  int ij, kl, ijkl;
  for (int i=0; i<nci; i++){
  for (int j=0; j<nci; j++){
  for (int k=0; k<nci; k++){
  for (int l=0; l<nci; l++){
      ij = INDEX(ndocc + i, ndocc + j);
      kl = INDEX(ndocc + k, ndocc + l);
      ijkl = INDEX(ij, kl);
      actTPDMp[i * nci3 + j * nci2 + k * nci + l] = fullTPDMp[ijkl];
  }}}}
  fullTPDM.reset();
  return actTPDM;
}
void CIWavefunction::set_tpdm()
{
  int npop = CalcInfo.num_ci_orbs + CalcInfo.num_drc_orbs;

  SharedVector symm_tpdm = get_tpdm(true);
  double* symm_tpdmp = symm_tpdm->pointer();

  SharedMatrix sm_tpdm(new Matrix("TPDM", npop * npop, npop * npop));
  double* sm_tpdmp = sm_tpdm->pointer()[0];

  int npop2 = npop * npop;
  int npop3 = npop * npop * npop;
  int ij, kl, ijkl;
  for (int i=0; i<npop; i++){
  for (int j=0; j<npop; j++){
  for (int k=0; k<npop; k++){
  for (int l=0; l<npop; l++){
      ij = INDEX(i,j);
      kl = INDEX(k,l);
      ijkl = INDEX(ij, kl);
      sm_tpdmp[i * npop3 + j * npop2 + k * npop + l] = symm_tpdmp[ijkl];
  }}}}
  symm_tpdm.reset();

  // Set numpy shape
  sm_tpdm->set_numpy_dims(4);
  int* shape = new int[4];
  shape[0] = npop; shape[1] = npop;
  shape[2] = npop; shape[3] = npop;
  sm_tpdm->set_numpy_shape(shape);

  TPDM_ = sm_tpdm;
}

// void CIWavefunction::finalize()
// {
//
// }


}}
