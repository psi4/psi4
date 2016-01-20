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

#include <libmints/mints.h>
#include <psi4-dec.h>
#include <libciomr/libciomr.h>
#include <psifiles.h>
#include <libqt/qt.h>
#include <libscf_solver/hf.h>

#include "globaldefs.h"
#include "ciwave.h"
#include "civect.h"
#include "structs.h"
#include "slaterd.h"

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
    Ca_ = reference_wavefunction_->Ca()->clone();
    Cb_ = Ca_; // We can only do RHF or ROHF reference wavefunctions.
    Da_ = reference_wavefunction_->Da()->clone(); // This will only be overwritten if form_opdm is called
    Db_ = reference_wavefunction_->Db()->clone();

    // Set information
    ints_init_ = false;
    df_ints_init_ = false;

    // Form strings
    outfile->Printf("\n   ==> Setting up CI strings <==\n\n");
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
    else if (orbitals == "ROT"){
      for (int h=0; h<nirrep_; h++){
        start[h] = CalcInfo_->frozen_docc[h];
        end[h] = nmopi_[h] - CalcInfo_->frozen_uocc[h];
      }
    }
    else if (orbitals == "ALL"){
      for (int h=0; h<nirrep_; h++){
        start[h] = 0;
        end[h] = nmopi_[h];
      }
    }
    else{
        throw PSIEXCEPTION("CIWave: Orbital subset is not defined, should be FZC, DRC, DOCC, ACT, RAS1, RAS2, RAS3, RAS4, POP, VIR, FZV, DRV, or ALL");
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

        if (opdm_map_.count(opdm_name.str()) == 0){
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
SharedMatrix CIWavefunction::get_tpdm(const std::string& spin, bool symmetrize)
{

  if (symmetrize){
    if (spin != "SUM")
        throw PSIEXCEPTION("CIWavefunction::get_tpdm: Symmetrize is only available for SUM spin type.");

    // Build
    int nact = CalcInfo_->num_ci_orbs;
    int nact2 = nact * nact;

    double** tpdm_nsp = tpdm_->pointer();
    SharedMatrix ret(new Matrix("MO-basis TPDM (symmetrized)", nact2, nact2));
    double** retp = ret->pointer();

    // Symmetrize
    for (int p=0, target=0; p<CalcInfo_->num_ci_orbs; p++) {
    for (int q=0; q<=p; q++) {
    for (int r=0; r<=p; r++) {

      int smax = (p == r) ? q+1 : r+1;
      for (int s=0; s<smax; s++) {

        // tpdm_nsp indices
       int pq = p * nact + q;
       int qp = q * nact + p;
       int rs = r * nact + s;
       int sr = s * nact + r;

       /* would be 0.25 but the formulae I used for the diag hessian
        * seem to define the TPDM with the 1/2 back outside */
       double value = 0.5 * (tpdm_nsp[pq][rs] + tpdm_nsp[qp][rs] +
                             tpdm_nsp[pq][sr] + tpdm_nsp[qp][sr]);

       // Write out 8 fold symmetry
       retp[pq][rs] = retp[qp][rs] =
       retp[pq][sr] = retp[qp][sr] =
       retp[rs][pq] = retp[rs][qp] =
       retp[sr][pq] = retp[sr][qp] = value;

      }
    }}}

    // Return
    return ret;
  }
  else {
    if (spin == "SUM")    return tpdm_;
    else if (spin == "AA") return tpdm_aa_;
    else if (spin == "AB") return tpdm_ab_;
    else if (spin == "BB") return tpdm_bb_;
    else throw PSIEXCEPTION("CIWavefunction::get_tpdm: Spin type must be AA, AB, BB, or SUM.");
  }
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

    // Cleanup up MCSCF integral objects
    if (Parameters_->mcscf){
        jk_.reset();
        if (MCSCF_Parameters_->mcscf_type == "DF"){
            dferi_.reset();
        }
        else{
            rot_space_.reset();
            act_space_.reset();
            ints_.reset();
        }
    }
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


SharedCIVector CIWavefunction::new_civector(int maxnvect, int filenum, bool use_disk,
                                            bool buf_init)
{
   SharedCIVector civect(new CIvect(Parameters_->icore, maxnvect, (int)use_disk,
                         filenum, CIblks_, CalcInfo_, Parameters_, H0block_, buf_init));
   return civect;

}
SharedMatrix CIWavefunction::hamiltonian(void)
{
    BIGINT size = CIblks_->vectlen;
    double h_size_gb = (double)(8 * size * size) / 1E9;
    if (h_size_gb > 1){
        outfile->Printf("CIWave::Requsted size of the hamiltonian is %lf!\n", h_size_gb);
        throw PSIEXCEPTION("CIWave::hamiltonian: Size is too large for explicit hamiltonian build");
    }

    SharedMatrix H(new Matrix("CI Hamiltonian", (int)size, (int)size));
    double** Hp = H->pointer();

    CIvect Cvec(1, 1, 0, 0, CIblks_, CalcInfo_, Parameters_, H0block_);
    SlaterDeterminant I, J;
    int Iarel, Ialist, Ibrel, Iblist;
    for (int ii=0; ii<size; ii++) {
        Cvec.det2strings(ii, &Ialist, &Iarel, &Iblist, &Ibrel);
        I.set(CalcInfo_->num_alp_expl,
             alplist_[Ialist][Iarel].occs, CalcInfo_->num_bet_expl,
             betlist_[Iblist][Ibrel].occs);
        Hp[ii][ii] = matrix_element(&I, &I) + CalcInfo_->edrc;

        /* introduce symmetry or other restrictions here */
        for (int jj=0; jj<ii; jj++) {
            Cvec.det2strings(jj, &Ialist, &Iarel, &Iblist, &Ibrel);
            J.set(CalcInfo_->num_alp_expl,
               alplist_[Ialist][Iarel].occs, CalcInfo_->num_bet_expl,
               betlist_[Iblist][Ibrel].occs);
            Hp[ii][jj] = Hp[jj][ii] = matrix_element(&I, &J);
        }
    }
    return H;
}
SharedMatrix CIWavefunction::orbital_ci_block(int fi, int fj)
{
 int filenum = Parameters_->d_filenum;
 CIvect Ivec(Parameters_->icore, 1, 1, filenum, CIblks_, CalcInfo_, Parameters_,
               H0block_, true);
 Ivec.init_io_files(true);
 Ivec.read(0, 0);

 // SharedVector oc_block(new Vector("OC Block", size));
 // double* oc_blockp = oc_block->pointer();
 int nci = CalcInfo_->num_ci_orbs;
 SharedMatrix ret(new Matrix("OPDM A Scratch", nci, nci));
 double** retp = ret->pointer();


  int Ia_idx, Ib_idx, Ja_idx, Jb_idx, Ja_ex, Jb_ex, Jbcnt, Jacnt;
  struct stringwr *Jb, *Ja;
  signed char *Jbsgn, *Jasgn;
  unsigned int *Jbridx, *Jaridx;
  double C1, C2, Ib_sgn, Ia_sgn;
  int i, j, oij, ndrc, *Jboij, *Jaoij;

  outfile->Printf("\n\nStarting my OPDM code\n");
  for (int Iblock=0; Iblock<Ivec.num_blocks_; Iblock++) {
    int Iac = Ivec.Ia_code_[Iblock];
    int Ibc = Ivec.Ib_code_[Iblock];
    int Inas = Ivec.Ia_size_[Iblock];
    int Inbs = Ivec.Ib_size_[Iblock];
    outfile->Printf("Iac = %d, Ibc = %d, Inas = %d, Inbs = %d\n", Iac, Ibc, Inas, Inbs);
    if (Inas==0 || Inbs==0) continue;
    for (int Jblock=0; Jblock<Ivec.num_blocks_; Jblock++) {
      int Jac = Ivec.Ia_code_[Jblock];
      int Jbc = Ivec.Ib_code_[Jblock];
      int Jnas = Ivec.Ia_size_[Jblock];
      int Jnbs = Ivec.Ib_size_[Jblock];
      if (!(s1_contrib_[Iblock][Jblock] || s2_contrib_[Iblock][Jblock])) continue;
      outfile->Printf("Jac = %d, Jbc = %d, Jnas = %d, Jnbs = %d\n", Jac, Jbc, Jnas, Jnbs);

//         opdm_block(alplist_, betlist_, scratch_ap, scratch_bp, Jvec.blocks_[Jblock],
//                    Ivec.blocks_[Iblock], Jac, Jbc, Jnas,
//                    Jnbs, Iac, Ibc, Inas, Inbs);
// void CIWavefunction::opdm_block(struct stringwr **alplist, struct stringwr **betlist,
//     double **onepdm_a, double **onepdm_b, double **CJ, double **CI, int Ja_list,
//     int Jb_list, int Jnas, int Jnbs, int Ia_list, int Ib_list,
//     int Inas, int Inbs)

      /* loop over Ia in Iac */
      if (Iac == Jac) {
        for (Ia_idx=0; Ia_idx<Inas; Ia_idx++) {
        for (Jb=betlist_[Jbc], Jb_idx=0; Jb_idx<Jnbs; Jb_idx++, Jb++) {
        C1 = Ivec.blocks_[Jblock][Ia_idx][Jb_idx];

        /* loop over excitations E^b_{ij} from |B(J_b)> */
        Jbcnt = Jb->cnt[Ibc];
        Jbridx = Jb->ridx[Ibc];
        Jbsgn = Jb->sgn[Ibc];
        Jboij = Jb->oij[Ibc];
        for (Jb_ex=0; Jb_ex < Jbcnt; Jb_ex++) {
          oij = *Jboij++;
          Ib_idx = *Jbridx++;
          Ib_sgn = (double) *Jbsgn++;
          C2 = Ivec.blocks_[Iblock][Ia_idx][Ib_idx];
                i = oij/CalcInfo_->num_ci_orbs;
                j = oij%CalcInfo_->num_ci_orbs;
          outfile->Printf("%d %d | %d | %d %d | %lf %lf %lf\n", Ia_idx, Jb_idx, Jb_ex, i, j, C1, C2, Ib_sgn);
          retp[i][j] += C1 * C2 * Ib_sgn;
        }
        outfile->Printf("---\n");
        }
        }
      }

      /* loop over Ib in Ibc */
      if (Ibc == Jbc) {
        for (Ib_idx=0; Ib_idx<Inbs; Ib_idx++) {
        for (Ja=alplist_[Jac], Ja_idx=0; Ja_idx<Jnas; Ja_idx++, Ja++) {
        C1 = Ivec.blocks_[Jblock][Ja_idx][Ib_idx];

        /* loop over excitations */
        Jacnt = Ja->cnt[Iac];
        Jaridx = Ja->ridx[Iac];
        Jasgn = Ja->sgn[Iac];
        Jaoij = Ja->oij[Iac];
        for (Ja_ex=0; Ja_ex < Jacnt; Ja_ex++) {
          oij = *Jaoij++;
          Ia_idx = *Jaridx++;
          Ia_sgn = (double) *Jasgn++;
          C2 = Ivec.blocks_[Iblock][Ia_idx][Ib_idx];
                i = oij/CalcInfo_->num_ci_orbs;
                j = oij%CalcInfo_->num_ci_orbs;
          retp[i][j] += C1 * C2 * Ia_sgn;
        }
        }
        }
      }

    } /* end loop over Jblock */
  } /* end loop over Iblock */
  return ret;

}

}} // End Psi and CIWavefunction spaces
