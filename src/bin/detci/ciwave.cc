#include <libmints/mints.h>
#include <psi4-dec.h>
#include <libiwl/iwl.h>
#include <libciomr/libciomr.h>
#include <psifiles.h>
#include <libqt/qt.h>

#include "globaldefs.h"
#include "ciwave.h"
#include "structs.h"


namespace psi { namespace detci {

CIWavefunction::CIWavefunction(boost::shared_ptr<Wavefunction> ref_wfn)
    : Wavefunction(Process::environment.options, ref_wfn->psio())
{
    set_reference_wavefunction(ref_wfn);
    common_init();
}

CIWavefunction::CIWavefunction(boost::shared_ptr<Wavefunction> ref_wfn,
                               Options &options)
    : Wavefunction(options, ref_wfn->psio())
{
    set_reference_wavefunction(ref_wfn);
    common_init();
}
CIWavefunction::~CIWavefunction()
{
}
void CIWavefunction::common_init()
{
    // Copy the wavefuntion then update
    copy(reference_wavefunction_);

    title();
    init_ioff();

    // Build and set structs
    sme_first_call_ = 1;
    MCSCF_Parameters_ = new mcscf_params();
    CIblks_ = new ci_blks();
    SigmaData_ = new sigma_data();
    CalcInfo_ = new calcinfo();
    Parameters_ = new params();
    H0block_ = new H_zero_block();

    // CI Params
    get_parameters(options_);     /* get running params (convergence, etc)    */
    get_mo_info();               /* read DOCC, SOCC, frozen, nmo, etc        */
    set_ras_parameters();             /* set fermi levels and the like            */

    print_parameters();
    print_ras_parameters();

    //MCSCF_Params
    get_mcscf_parameters();

    // Wavefunction frozen nomenclature is equivalent to dropped in detci.
    // In detci frozen means doubly occupied, but no orbital rotations.
    nalpha_     = CalcInfo_->num_alp; // Total number of alpha electrons, including core
    nbeta_      = CalcInfo_->num_bet;
    nfrzc_      = CalcInfo_->num_drc_orbs;

    // Per irrep data, this is approximate and should typically not be used
    // Data should come from get_dimension(space)
    for(int h = 0; h < nirrep_; ++h){
        frzcpi_[h] = CalcInfo_->dropped_docc[h];
        doccpi_[h] = CalcInfo_->docc[h];
        soccpi_[h] = CalcInfo_->socc[h];
        frzvpi_[h] = CalcInfo_->dropped_uocc[h];

        nmopi_[h]  = CalcInfo_->orbs_per_irr[h];
        nsopi_[h]  = CalcInfo_->so_per_irr[h];
    }

    // Set relevant matrices
    //H_ = reference_wavefunction_->H();
    Ca_ = reference_wavefunction_->Ca()->clone();
    Cb_ = Ca_; // We can only do RHF or ROHF reference wavefunctions.
    Da_ = reference_wavefunction_->Da()->clone(); // This will only be overwritten if form_opdm is called
    Db_ = reference_wavefunction_->Db()->clone();

    // Set information
    ints_init_ = false;
    df_ints_init_ = false;

    // Form strings
    form_strings();

    // Form Bendazzoli OV arrays
    if (Parameters_->bendazzoli) form_ov();

    name_ = "CIWavefunction";
}
double CIWavefunction::compute_energy()
{

   if (Parameters_->istop) {      /* Print size of space, other stuff, only   */
     cleanup();
     Process::environment.globals["CURRENT ENERGY"] = 0.0;
     Process::environment.globals["CURRENT CORRELATION ENERGY"] = 0.0;
     Process::environment.globals["CI TOTAL ENERGY"] = 0.0;
     Process::environment.globals["CI CORRELATION ENERGY"] = 0.0;

     return Success;
   }

   // MCSCF is special, we let it handle a lot of its own issues
   if (Parameters_->mcscf){
     compute_mcscf();
   }
   else{
     // Transform and set ci integrals
     transform_ci_integrals();

     if (Parameters_->mpn){
       compute_mpn();
       }
     else if (Parameters_->cc)
       compute_cc();
     else
       diag_h();
   }

   // Finished CI, setting wavefunction parameters
   if(!Parameters_->zaptn & Parameters_->opdm){
     form_opdm();
   }

   if (Parameters_->dipmom) opdm_properties();
   if (Parameters_->opdm_diag) ci_nat_orbs();
   if (Parameters_->tpdm) form_tpdm();
   if (Parameters_->print_lvl > 0){
     outfile->Printf("\t\t \"A good bug is a dead bug\" \n\n");
     outfile->Printf("\t\t\t - Starship Troopers\n\n");
     outfile->Printf("\t\t \"I didn't write FORTRAN.  That's the problem.\"\n\n");
     outfile->Printf("\t\t\t - Edward Valeev\n\n");
   }

   cleanup();

   return Process::environment.globals["CURRENT ENERGY"];
}

void CIWavefunction::orbital_locations(const std::string& orbitals, int* start, int* end){
    if (orbitals == "FZC"){
      for (int h=0; h<nirrep_; h++){
        start[h] = 0;
        end[h] = CalcInfo_->frozen_docc[h];
      }
    }
    else if (orbitals == "DOCC"){
      for (int h=0; h<nirrep_; h++){
        start[h] = CalcInfo_->frozen_docc[h];
        end[h] = CalcInfo_->dropped_docc[h];
      }
    }
    else if (orbitals == "DRC"){
      for (int h=0; h<nirrep_; h++){
        start[h] = 0;
        end[h] = CalcInfo_->dropped_docc[h];
      }
    }
    else if (orbitals == "ACT"){
      for (int h=0; h<nirrep_; h++){
        start[h] = CalcInfo_->dropped_docc[h];
        end[h] = nmopi_[h] - CalcInfo_->dropped_uocc[h];
      }
    }
    else if (orbitals == "RAS1"){
      for (int h=0; h<nirrep_; h++){
        start[h] = CalcInfo_->dropped_docc[h];
        end[h] = start[h] + CalcInfo_->ras_opi[0][h];
      }
    }
    else if (orbitals == "RAS2"){
      for (int h=0; h<nirrep_; h++){
        start[h] = CalcInfo_->dropped_docc[h] + CalcInfo_->ras_opi[0][h];
        end[h] = start[h] + CalcInfo_->ras_opi[1][h];
      }
    }
    else if (orbitals == "RAS3"){
      for (int h=0; h<nirrep_; h++){
        start[h] = CalcInfo_->dropped_docc[h] + CalcInfo_->ras_opi[0][h] + CalcInfo_->ras_opi[1][h];
        end[h] = start[h] + CalcInfo_->ras_opi[2][h];
      }
    }
    else if (orbitals == "RAS4"){
      for (int h=0; h<nirrep_; h++){
        start[h] = nmopi_[h] - CalcInfo_->dropped_uocc[h] - CalcInfo_->ras_opi[3][h];
        end[h] = nmopi_[h] - CalcInfo_->dropped_uocc[h];
      }
    }
    else if (orbitals == "POP"){
      for (int h=0; h<nirrep_; h++){
        start[h] = 0;
        end[h] = nmopi_[h] - CalcInfo_->dropped_uocc[h];
      }
    }
    else if (orbitals == "DRV"){
      for (int h=0; h<nirrep_; h++){
        start[h] = nmopi_[h] - CalcInfo_->dropped_uocc[h];
        end[h] = nmopi_[h];
      }
    }
    else if (orbitals == "VIR"){
      for (int h=0; h<nirrep_; h++){
        start[h] = nmopi_[h] - CalcInfo_->dropped_uocc[h];
        end[h] = nmopi_[h] - CalcInfo_->frozen_uocc[h];
      }
    }
    else if (orbitals == "FZV"){
      for (int h=0; h<nirrep_; h++){
        start[h] = nmopi_[h] - CalcInfo_->frozen_uocc[h];
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
        start[h] = CalcInfo_->frozen_docc[h];
        end[h] = nmopi_[h] - CalcInfo_->frozen_uocc[h];
      }
    }
    else{
        throw PSIEXCEPTION("CIWAVE: Orbital subset is not defined, should be FZC, DRC, DOCC, ACT, RAS1, RAS2, RAS3, RAS4, POP, VIR, FZV, DRV, or ALL");
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

Dimension CIWavefunction::get_dimension(const std::string& orbital_name)
{
    /// Figure out orbital positions
    int* start = new int[nirrep_];
    int* end = new int[nirrep_];
    orbital_locations(orbital_name, start, end);

    Dimension dim = Dimension(nirrep_);

    for (int h=0; h<nirrep_; h++){
      dim[h] = end[h] - start[h];
    }

    delete[] start;
    delete[] end;
    return dim;
}

SharedMatrix CIWavefunction::get_opdm(int Iroot, int Jroot, const std::string& spin, bool full_space)
{

    double inact_value = (spin == "SUM") ? 2.0 : 1.0;
    SharedMatrix opdm;

    if ((Iroot == -1) && (Jroot == -1)){
        if (spin == "SUM") opdm = opdm_;
        else if (spin == "A") opdm = opdm_a_;
        else if (spin == "B") opdm = opdm_b_;
        else throw PSIEXCEPTION("CIWavefunction::get_opdm: Spin type must be A, B, or SUM.");
    }
    else {
        if (Jroot == -1) Jroot = Iroot;

        std::stringstream opdm_name;
        if (spin == "SUM") opdm_name << "MO-basis OPDM <" << Iroot+1 << "| Etu |" << Jroot+1 << ">";
        else if (spin == "A") opdm_name << "MO-basis Alpha OPDM <" << Iroot+1 << "| Etu |" << Jroot+1 << ">";
        else if (spin == "B") opdm_name << "MO-basis Beta OPDM <" << Iroot+1 << "| Etu |" << Jroot+1 << ">";
        else throw PSIEXCEPTION("CIWavefunction::get_opdm: Spin type must be A, B, or SUM.");

        if (opdm_map_.find(opdm_name.str()) != opdm_map_.end()){
            throw PSIEXCEPTION("CIWavefunction::get_opdm: Requested OPDM was not formed!\n");
        }

        opdm = opdm_map_[opdm_name.str()];
    }

    if (full_space) {
        return opdm_add_inactive(opdm, inact_value, true);
    }
    else {
        return opdm;
    }
}
SharedVector CIWavefunction::get_tpdm(bool symmetrize, const std::string& tpdm_type)
{

  int sqnbf, ntri;
  int *ioff_lt, i;                    /* offsets for left (or right) indices */
  int p,q,r,s,smax,pq,qp,rs,sr,pqrs,qprs,pqsr,qpsr,target;
  struct iwlbuf TBuff;

  int tpdm_file;
  if (tpdm_type == "SUM") tpdm_file = Parameters_->tpdm_file;
  else if (tpdm_type == "AA") tpdm_file = PSIF_MO_AA_TPDM;
  else if (tpdm_type == "BB") tpdm_file = PSIF_MO_BB_TPDM;
  else if (tpdm_type == "AB") tpdm_file = PSIF_MO_AB_TPDM;
  else throw PSIEXCEPTION("CIWAVE: TPDM can only be SUM, AA, AB, or BB");

  iwl_buf_init(&TBuff, tpdm_file, 0.0, 1, 1);
  int npop = CalcInfo_->num_ci_orbs + CalcInfo_->num_drc_orbs;

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

 iwl_buf_rd_all(&TBuff, tpdmp, ioff_lt, ioff_lt, 1, ioff_,
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
SharedMatrix CIWavefunction::get_active_tpdm(const std::string& tpdm_type)
{
  int nci   = CalcInfo_->num_ci_orbs;
  int ndocc = CalcInfo_->num_drc_orbs;
  int npop  = CalcInfo_->num_ci_orbs + CalcInfo_->num_drc_orbs;
  SharedMatrix actTPDM(new Matrix("TPDM", nci*nci, nci*nci));
  double* actTPDMp = actTPDM->pointer()[0];

  SharedVector fullTPDM = get_tpdm(true, tpdm_type);
  double* fullTPDMp = fullTPDM->pointer();

  int nci2 = nci * nci;
  int nci3 = nci * nci * nci;
  for (int i=0; i<nci; i++){
  for (int j=0; j<nci; j++){
  for (int k=0; k<nci; k++){
  for (int l=0; l<nci; l++){
      int r_i = CalcInfo_->act_reorder[i];
      int r_j = CalcInfo_->act_reorder[j];
      int r_k = CalcInfo_->act_reorder[k];
      int r_l = CalcInfo_->act_reorder[l];
      int r_ij = INDEX(ndocc + r_i, ndocc + r_j);
      int r_kl = INDEX(ndocc + r_k, ndocc + r_l);
      int r_ijkl = INDEX(r_ij, r_kl);
      actTPDMp[i * nci3 + j * nci2 + k * nci + l] = fullTPDMp[r_ijkl];
  }}}}

  fullTPDM.reset();

  // Set numpy shape
  actTPDM->set_numpy_dims(4);
  int* shape = new int[4];
  shape[0] = nci; shape[1] = nci;
  shape[2] = nci; shape[3] = nci;
  actTPDM->set_numpy_shape(shape);

  return actTPDM;
}


void CIWavefunction::form_tpdm(void)
{
  tpdm(alplist_, betlist_, Parameters_->num_roots,
       Parameters_->num_d_tmp_units, Parameters_->first_d_tmp_unit,
       Parameters_->num_roots,
       Parameters_->num_d_tmp_units, Parameters_->first_d_tmp_unit,
       Parameters_->tpdm_file, Parameters_->tpdm_write, Parameters_->tpdm_print);
}

/*
** cleanup(): Free any allocated memory that wasn't already freed elsewhere
*/
void CIWavefunction::cleanup(void)
{
    delete[] ioff_;
    sigma_free();

    // Free Bendazzoli OV arrays
    //if (Parameters_->bendazzoli) free(OV);

    // CalcInfo free
    free(CalcInfo_->onel_ints);
    free(CalcInfo_->twoel_ints);
    free(CalcInfo_->maxK);
    free_matrix(CalcInfo_->gmat, CalcInfo_->num_ci_orbs);

    // DGAS main areas to track size of
    // Free strings and graphs
    //delete SigmaData_;

    //delete MCSCF_Parameters_;
    //delete CIblks_;
    //delete CalcInfo_;
    //delete Parameters_;
    //delete H0block_;
}

/*
** title(): Function prints a program identification
*/
void CIWavefunction::title(void)
{
  // if (Parameters_->print_lvl) {
   outfile->Printf("\n");
   outfile->Printf("         ---------------------------------------------------------\n");
   outfile->Printf("                                 D E T C I  \n");
   outfile->Printf("\n");
   outfile->Printf("                             C. David Sherrill\n") ;
   outfile->Printf("                             Matt L. Leininger\n") ;
   outfile->Printf("                               18 June 1999\n") ;
   outfile->Printf("         ---------------------------------------------------------\n");
   outfile->Printf("\n");
}


/*
** init_ioff(): Set up the ioff array for quick indexing
**
*/
void CIWavefunction::init_ioff(void)
{
   ioff_ = new int[65025];
   ioff_[0] = 0;
   for (int i = 1; i < 65025; i++) {
      ioff_[i] = ioff_[i-1] + i;
    }
}


// void CIWavefunction::finalize()
// {
//
// }


}}
