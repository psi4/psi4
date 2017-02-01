/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
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
 * @END LICENSE
 */

 #include "psi4/pragma.h"
 PRAGMA_WARNING_PUSH
 PRAGMA_WARNING_IGNORE_DEPRECATED_DECLARATIONS
 #include <memory>
 PRAGMA_WARNING_POP
#include "psi4/psi4-dec.h"
#include "psi4/libciomr/libciomr.h"
#include "psi4/psifiles.h"
#include "psi4/libqt/qt.h"
#include "psi4/libscf_solver/hf.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libfock/soscf.h"
#include "psi4/detci/globaldefs.h"
#include "psi4/detci/ciwave.h"
#include "psi4/detci/civect.h"
#include "psi4/detci/structs.h"
#include "psi4/detci/slaterd.h"

namespace psi { namespace detci {

CIWavefunction::CIWavefunction(std::shared_ptr<Wavefunction> ref_wfn)
    : Wavefunction(Process::environment.options) {
    // Copy the wavefuntion then update
    shallow_copy(ref_wfn);
    set_reference_wavefunction(ref_wfn);
    common_init();
}

CIWavefunction::CIWavefunction(std::shared_ptr<Wavefunction> ref_wfn,
                               Options& options) : Wavefunction(options) {

    // Copy the wavefuntion then update
    shallow_copy(ref_wfn);
    set_reference_wavefunction(ref_wfn);
    common_init();
}

CIWavefunction::~CIWavefunction() {
    cleanup_ci();
    cleanup_dpd();
}

void CIWavefunction::common_init() {
    title((options_.get_str("WFN") == "CASSCF") || (options_.get_str("WFN") == "RASSCF"));

    // Build and set structs
    sme_first_call_ = 1;
    CIblks_ = new ci_blks();
    SigmaData_ = new sigma_data();
    CalcInfo_ = new calcinfo();
    Parameters_ = new params();
    H0block_ = new H_zero_block();

    // CI Params
    get_parameters(options_); /* get running params (convergence, etc)    */
    get_mo_info();            /* read DOCC, SOCC, frozen, nmo, etc        */
    set_ras_parameters();     /* set fermi levels and the like            */

    // Print out information
    print_parameters();
    print_ras_parameters();

    // Wavefunction frozen nomenclature is equivalent to dropped in detci.
    // In detci frozen means doubly occupied, but no orbital rotations.
    nalpha_ = CalcInfo_->num_alp;  // Total number of alpha electrons, including core
    nbeta_ = CalcInfo_->num_bet;
    nfrzc_ = CalcInfo_->num_drc_orbs;

    // Per irrep data, this is approximate and should typically not be used
    // Data should come from get_dimension(space)
    for (int h = 0; h < nirrep_; ++h) {
        frzcpi_[h] = CalcInfo_->dropped_docc[h];
        doccpi_[h] = CalcInfo_->docc[h];
        soccpi_[h] = CalcInfo_->socc[h];
        frzvpi_[h] = CalcInfo_->dropped_uocc[h];
    }

    // Set relevant matrices
    Ca_ = reference_wavefunction_->Ca()->clone();
    Cb_ = Ca_;  // We can only do RHF or ROHF reference wavefunctions.
    Da_ = reference_wavefunction_->Da()->clone();
    Db_ = reference_wavefunction_->Db()->clone();

    // Set information
    ints_init_ = false;
    df_ints_init_ = false;
    mcscf_object_init_ = false;
    cleaned_up_ci_ = false;
    fzc_fock_computed_ = false;

    // Form strings
    outfile->Printf("\n   ==> Setting up CI strings <==\n\n");
    form_strings();

    // Form Bendazzoli OV arrays
    if (Parameters_->bendazzoli) form_ov();

    name_ = "CIWavefunction";

    // Init H0 block
    H0block_init(CIblks_->vectlen);

    if (CIblks_->vectlen < 2){
        throw PSIEXCEPTION("CIWavefunction: Must have more than one determinant!");
    }
}
size_t CIWavefunction::ndet() { return (size_t)CIblks_->vectlen; }

double CIWavefunction::compute_energy() {
    if (Parameters_->istop) { /* Print size of space, other stuff, only   */
        cleanup_ci();
        cleanup_dpd();
        Process::environment.globals["CURRENT ENERGY"] = 0.0;
        Process::environment.globals["CURRENT CORRELATION ENERGY"] = 0.0;
        Process::environment.globals["CI TOTAL ENERGY"] = 0.0;
        Process::environment.globals["CI CORRELATION ENERGY"] = 0.0;

        return 0.0;
    }

    // Transform and set ci integrals
    transform_ci_integrals();

    if (Parameters_->mpn) {
        compute_mpn();
    } else if (Parameters_->cc) {
        compute_cc();
    } else {
        diag_h();
    }

    // Finished CI, setting wavefunction parameters
    if (!Parameters_->zaptn &&
        (Parameters_->opdm || Parameters_->transdens || Parameters_->opdm_diag)) {
        form_opdm();
    }

    // if (Parameters_->opdm_diag) ci_nat_orbs();
    if (Parameters_->tpdm) form_tpdm();

    // cleanup_ci();
    // cleanup_dpd();

    return Process::environment.globals["CURRENT ENERGY"];
}

void CIWavefunction::orbital_locations(const std::string& orbitals, int* start,
                                       int* end) {
    if (orbitals == "FZC") {
        for (int h = 0; h < nirrep_; h++) {
            start[h] = 0;
            end[h] = CalcInfo_->frozen_docc[h];
        }
    } else if (orbitals == "DOCC") {
        for (int h = 0; h < nirrep_; h++) {
            start[h] = CalcInfo_->frozen_docc[h];
            end[h] = CalcInfo_->dropped_docc[h];
        }
    } else if (orbitals == "DRC") {
        for (int h = 0; h < nirrep_; h++) {
            start[h] = 0;
            end[h] = CalcInfo_->dropped_docc[h];
        }
    } else if (orbitals == "ACT") {
        for (int h = 0; h < nirrep_; h++) {
            start[h] = CalcInfo_->dropped_docc[h];
            end[h] = nmopi_[h] - CalcInfo_->dropped_uocc[h];
        }
    } else if (orbitals == "RAS1") {
        for (int h = 0; h < nirrep_; h++) {
            start[h] = CalcInfo_->dropped_docc[h];
            end[h] = start[h] + CalcInfo_->ras_opi[0][h];
        }
    } else if (orbitals == "RAS2") {
        for (int h = 0; h < nirrep_; h++) {
            start[h] = CalcInfo_->dropped_docc[h] + CalcInfo_->ras_opi[0][h];
            end[h] = start[h] + CalcInfo_->ras_opi[1][h];
        }
    } else if (orbitals == "RAS3") {
        for (int h = 0; h < nirrep_; h++) {
            start[h] = CalcInfo_->dropped_docc[h] + CalcInfo_->ras_opi[0][h] + CalcInfo_->ras_opi[1][h];
            end[h] = start[h] + CalcInfo_->ras_opi[2][h];
        }
    } else if (orbitals == "RAS4") {
        for (int h = 0; h < nirrep_; h++) {
            start[h] = nmopi_[h] - CalcInfo_->dropped_uocc[h] - CalcInfo_->ras_opi[3][h];
            end[h] = nmopi_[h] - CalcInfo_->dropped_uocc[h];
        }
    } else if (orbitals == "POP") {
        for (int h = 0; h < nirrep_; h++) {
            start[h] = 0;
            end[h] = nmopi_[h] - CalcInfo_->dropped_uocc[h];
        }
    } else if (orbitals == "DRV") {
        for (int h = 0; h < nirrep_; h++) {
            start[h] = nmopi_[h] - CalcInfo_->dropped_uocc[h];
            end[h] = nmopi_[h];
        }
    } else if (orbitals == "VIR") {
        for (int h = 0; h < nirrep_; h++) {
            start[h] = nmopi_[h] - CalcInfo_->dropped_uocc[h];
            end[h] = nmopi_[h] - CalcInfo_->frozen_uocc[h];
        }
    } else if (orbitals == "FZV") {
        for (int h = 0; h < nirrep_; h++) {
            start[h] = nmopi_[h] - CalcInfo_->frozen_uocc[h];
            end[h] = nmopi_[h];
        }
    } else if (orbitals == "ROT") {
        for (int h = 0; h < nirrep_; h++) {
            start[h] = CalcInfo_->frozen_docc[h];
            end[h] = nmopi_[h] - CalcInfo_->frozen_uocc[h];
        }
    } else if (orbitals == "OA") {
        for (int h = 0; h < nirrep_; h++) {
            start[h] = CalcInfo_->frozen_docc[h];
            end[h] = nmopi_[h] - CalcInfo_->dropped_uocc[h];
        }
    } else if (orbitals == "AV") {
        for (int h = 0; h < nirrep_; h++) {
            start[h] = CalcInfo_->dropped_docc[h];
            end[h] = nmopi_[h] - CalcInfo_->frozen_uocc[h];
        }
    } else if (orbitals == "ALL") {
        for (int h = 0; h < nirrep_; h++) {
            start[h] = 0;
            end[h] = nmopi_[h];
        }
    } else {
        throw PSIEXCEPTION(
            "CIWave: Orbital subset is not defined, should be FZC, DRC, DOCC, "
            "ACT, RAS1, RAS2, RAS3, RAS4, POP, VIR, FZV, DRV, OA, AV, ROT, or ALL");
    }
}

SharedMatrix CIWavefunction::get_orbitals(const std::string& orbital_name) {
    /// Figure out orbital positions
    int* start = new int[nirrep_];
    int* end = new int[nirrep_];

    orbital_locations(orbital_name, start, end);

    int* spread = new int[nirrep_];
    for (int h = 0; h < nirrep_; h++) {
        spread[h] = end[h] - start[h];
    }

    /// Fill desired orbitals
    SharedMatrix retC(new Matrix("C " + orbital_name, nirrep_, nsopi_, spread));
    for (int h = 0; h < nirrep_; h++) {
        for (int i = start[h], pos = 0; i < end[h]; i++, pos++) {
            C_DCOPY(nsopi_[h], &Ca_->pointer(h)[0][i], nmopi_[h],
                    &retC->pointer(h)[0][pos], spread[h]);
        }
    }

    /// Cleanup
    delete[] start;
    delete[] end;
    delete[] spread;

    return retC;
}

void CIWavefunction::set_orbitals(const std::string& orbital_name,
                                  SharedMatrix orbitals) {
    /// Figure out orbital positions
    int* start = new int[nirrep_];
    int* end = new int[nirrep_];

    orbital_locations(orbital_name, start, end);

    int* spread = new int[nirrep_];
    for (int h = 0; h < nirrep_; h++) {
        spread[h] = end[h] - start[h];
    }

    /// Fill desired orbitals
    for (int h = 0; h < nirrep_; h++) {
        for (int i = start[h], pos = 0; i < end[h]; i++, pos++) {
            C_DCOPY(nsopi_[h], &orbitals->pointer(h)[0][pos], spread[h],
                    &Ca_->pointer(h)[0][i], nmopi_[h]);
        }
    }

    /// Cleanup
    delete[] start;
    delete[] end;
    delete[] spread;
}

Dimension CIWavefunction::get_dimension(const std::string& orbital_name) {
    /// Figure out orbital positions
    int* start = new int[nirrep_];
    int* end = new int[nirrep_];
    orbital_locations(orbital_name, start, end);

    Dimension dim = Dimension(nirrep_);

    for (int h = 0; h < nirrep_; h++) {
        dim[h] = end[h] - start[h];
    }

    delete[] start;
    delete[] end;
    return dim;
}

SharedMatrix CIWavefunction::get_opdm(int Iroot, int Jroot,
                                      const std::string& spin,
                                      bool full_space) {
    if (!opdm_called_) {
        throw PSIEXCEPTION("CIWavefunction::get_opdm: OPDM was not formed!");
    }
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
        if (spin == "SUM") opdm_name << "MO-basis OPDM <" << Iroot << "| Etu |" << Jroot << ">";
        else if (spin == "A") opdm_name << "MO-basis Alpha OPDM <" << Iroot << "| Etu |" << Jroot << ">";
        else if (spin == "B") opdm_name << "MO-basis Beta OPDM <" << Iroot << "| Etu |" << Jroot << ">";
        else throw PSIEXCEPTION("CIWavefunction::get_opdm: Spin type must be A, B, or SUM.");

        if (opdm_map_.count(opdm_name.str()) == 0){
            throw PSIEXCEPTION("CIWavefunction::get_opdm: Requested OPDM was not formed!\n");
        }

        opdm = opdm_map_[opdm_name.str()];
    }

    if (Iroot != Jroot){ // Transition densities
        inact_value = 0.0;
    }

    if (full_space) {
        return opdm_add_inactive(opdm, inact_value, true);
    }
    else {
        return opdm;
    }
}
void CIWavefunction::convergence_death(){
    if (Parameters_->die_if_not_converged){
        throw PSIEXCEPTION("CIWavefunction: Iterations did not converge!");
    }

}
SharedMatrix CIWavefunction::get_tpdm(const std::string& spin,
                                      bool symmetrize) {
    if (!tpdm_called_) {
        throw PSIEXCEPTION("CIWavefunction::get_opdm: OPDM was not formed!");
    }

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
      for (int p=0; p<nact; p++) {
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

      // Add numpy shape
      std::vector<int> nshape{nact, nact, nact, nact};
      ret->set_numpy_shape(nshape);

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
void CIWavefunction::cleanup_ci(void) {

    // Make sure we dont double clean
    if (!cleaned_up_ci_){

        // Free Bendazzoli OV arrays
        // if (Parameters_->bendazzoli) free(OV);

        // DGAS main areas to track size of
        // Free strings and graphs

        // Free objects built in common_init
        if (CalcInfo_->sigma_initialized) sigma_free();
        delete SigmaData_;

        free_int_matrix(CIblks_->decode);
        free(CIblks_->first_iablk);
        free(CIblks_->last_iablk);
        delete CIblks_;

        //delete Parameters_;

        // Free H0block
        H0block_free();
        delete H0block_;

        // CalcInfo free
        free_int_matrix(CalcInfo_->ras_opi);
        for (int i = 0; i < 4; i++) {
            free_int_matrix(CalcInfo_->ras_orbs[i]);
        };

        cleaned_up_ci_ = true;
    }

}
void CIWavefunction::cleanup_dpd(void) {

    if (ints_init_) {
        ints_.reset();
        ints_init_ = false;
    }
    if (mcscf_object_init_){
        somcscf_.reset();
        mcscf_object_init_ = false;
    }
}
void CIWavefunction::set_ci_guess(std::string guess) {
    if (guess == "UNIT") {
        Parameters_->guess_vector = PARM_GUESS_VEC_UNIT;
    } else if (guess == "H0_BLOCK") {
        Parameters_->guess_vector = PARM_GUESS_VEC_H0_BLOCK;
    } else if (guess == "DFILE") {
        Parameters_->guess_vector = PARM_GUESS_VEC_DFILE;
    } else {
        throw PSIEXCEPTION(
            "CIWavefunction::set_ci_guess: Guess can only be UNIT, H0_BLOCK, or DFILE");
    }
}
void CIWavefunction::title(bool is_mcscf) {
    if (is_mcscf) {
        outfile->Printf("\n");
        outfile->Printf("         ---------------------------------------------------------\n");
        outfile->Printf("                Multi-Configurational Self-Consistent Field\n");
        outfile->Printf("                            (a 'D E T C I' module)\n");
        outfile->Printf("\n");
        outfile->Printf("                 Daniel G. A. Smith, C. David Sherrill, and\n");
        outfile->Printf("                              Matt L. Leininger\n");
        outfile->Printf("         ---------------------------------------------------------\n");
        outfile->Printf("\n");
    } else {
        outfile->Printf("\n");
        outfile->Printf("         ---------------------------------------------------------\n");
        outfile->Printf("                          Configuration Interaction\n");
        outfile->Printf("                            (a 'D E T C I' module)\n");
        outfile->Printf("\n");
        outfile->Printf("                 C. David Sherrill, Daniel G. A. Smith, and\n");
        outfile->Printf("                              Matt L. Leininger\n");
        outfile->Printf("         ---------------------------------------------------------\n");
        outfile->Printf("\n");
    }
}

SharedCIVector CIWavefunction::new_civector(int maxnvect, int filenum,
                                            bool use_disk, bool buf_init) {
    SharedCIVector civect(new CIvect(Parameters_->icore, maxnvect,
                                     (int)use_disk, filenum, CIblks_, CalcInfo_,
                                     Parameters_, H0block_, buf_init));
    return civect;
}
SharedCIVector CIWavefunction::D_vector(){
    SharedCIVector civect = new_civector(Parameters_->num_roots, Parameters_->d_filenum, true, true);
    return civect;
}
SharedCIVector CIWavefunction::Hd_vector(int hd_type) {
    hd_type = (hd_type == -1) ? Parameters_->hd_ave : hd_type;
    SharedCIVector Hd = new_civector(1, Parameters_->hd_filenum, true, true);
    Hd->init_io_files(false);  // False for do not open old
    Hd->diag_mat_els(alplist_, betlist_, CalcInfo_->onel_ints->pointer(),
                     CalcInfo_->twoel_ints->pointer(), CalcInfo_->edrc,
                     CalcInfo_->num_alp_expl, CalcInfo_->num_bet_expl,
                     CalcInfo_->num_ci_orbs, hd_type);
    Hd->write(0, 0);
    return Hd;
}
SharedMatrix CIWavefunction::hamiltonian(size_t hsize) {
    BIGINT size = (hsize) ? (BIGINT)hsize : CIblks_->vectlen;
    double h_size = (double)(8 * size * size);
    if (h_size > (Process::environment.get_memory() * 0.4)) {
        outfile->Printf("CIWave::Requsted size of the hamiltonian is %lf!\n", h_size / 1E9);
        throw PSIEXCEPTION("CIWave::hamiltonian: Size is too large for"
                           "explicit hamiltonian build");
    }

    SharedMatrix H(new Matrix("CI Hamiltonian", (size_t)size, (size_t)size));
    double** Hp = H->pointer();

    CIvect Cvec(1, 1, 0, 0, CIblks_, CalcInfo_, Parameters_, H0block_);

    SlaterDeterminant I, J;
    int Iarel, Ialist, Ibrel, Iblist;
    for (size_t ii = 0; ii < size; ii++) {
        Cvec.det2strings(ii, &Ialist, &Iarel, &Iblist, &Ibrel);
        I.set(CalcInfo_->num_alp_expl, alplist_[Ialist][Iarel].occs,
              CalcInfo_->num_bet_expl, betlist_[Iblist][Ibrel].occs);
        Hp[ii][ii] = matrix_element(&I, &I) + CalcInfo_->edrc;

        /* introduce symmetry or other restrictions here */
        for (size_t jj = 0; jj < ii; jj++) {
            Cvec.det2strings(jj, &Ialist, &Iarel, &Iblist, &Ibrel);
            J.set(CalcInfo_->num_alp_expl, alplist_[Ialist][Iarel].occs,
                  CalcInfo_->num_bet_expl, betlist_[Iblist][Ibrel].occs);
            Hp[ii][jj] = Hp[jj][ii] = matrix_element(&I, &J);
        }
    }
    return H;
    /* construct and print one block at a time for debugging */
    /*
    int ii2, jj2, blk, blk2, det1, det2;
    double **Hpart;

    for (blk = 0; blk < CIblks_->num_blocks; blk++) {
      for (blk2 = 0; blk2 < CIblks_->num_blocks; blk2++) {
        Hpart = init_matrix(CIblks_->Ia_size[blk]*CIblks_->Ib_size[blk],
                            CIblks_->Ia_size[blk2]*CIblks_->Ib_size[blk2]);
        for (ii=0,det1=0; ii<CIblks_->Ia_size[blk]; ii++) {
          for (jj=0; jj<CIblks_->Ib_size[blk]; jj++, det1++) {
            I.set(CalcInfo_->num_alp_expl,alplist[CIblks_->Ia_code[blk]][ii].occs,
                 CalcInfo_->num_bet_expl,betlist[CIblks_->Ib_code[blk]][jj].occs);
            for (ii2=0,det2=0; ii2<CIblks_->Ia_size[blk2]; ii2++) {
              for (jj2=0; jj2<CIblks_->Ib_size[blk2]; jj2++,det2++) {
                J.set(CalcInfo_->num_alp_expl,
                      alplist[CIblks_->Ia_code[blk2]][ii2].occs,
                      CalcInfo_->num_bet_expl,
                      betlist[CIblks_->Ib_code[blk2]][jj2].occs);
                Hpart[det1][det2] = matrix_element(&I,&J);
              }
            }
          }
        }
        if (print_ > 4 && size < 200) {
          outfile->Printf( "\nBlock %d %d of ", blk, blk2);
          outfile->Printf( "Hamiltonian matrix:\n");
          print_mat(Hpart, CIblks_->Ia_size[blk]*CIblks_->Ib_size[blk],
                           CIblks_->Ia_size[blk2]*CIblks_->Ib_size[blk2],
                    outfile);
        }
        free_matrix(Hpart, CIblks_->Ia_size[blk]*CIblks_->Ib_size[blk]);
      }
    }
    */
    /* end block-at-a-time stuff */

}

void CIWavefunction::init_mcscf_object(){

    if (Parameters_->mcscf_type == "DF") {
        if (!df_ints_init_) setup_dfmcscf_ints();
        somcscf_ = std::shared_ptr<SOMCSCF>(new DFSOMCSCF(jk_, dferi_, AO2SO_, H_));
    } else if (Parameters_->mcscf_type == "AO" ){
        if (!ints_init_) setup_mcscf_ints_ao();
        somcscf_ = std::shared_ptr<SOMCSCF>(new IncoreSOMCSCF(jk_, AO2SO_, H_));
    }
    else {
        if (!ints_init_) setup_mcscf_ints();
        somcscf_ = std::shared_ptr<SOMCSCF>(new DiskSOMCSCF(jk_, ints_, AO2SO_, H_));
    }

    // We assume some kind of ras here.
    if (Parameters_->wfn != "CASSCF") {
        std::vector<Dimension> ras_spaces;

        // We only have four spaces currently
        for (int nras = 0; nras < 4; nras++) {
            Dimension rasdim = Dimension(nirrep_, "RAS" + std::to_string(nras));
            for (int h = 0; h < nirrep_; h++) {
                rasdim[h] = CalcInfo_->ras_opi[nras][h];
            }
            ras_spaces.push_back(rasdim);
        }
        somcscf_->set_ras(ras_spaces);
    }

    // Set fzc energy
    SharedMatrix Cfzc = get_orbitals("FZC");
    somcscf_->set_frozen_orbitals(Cfzc);
    mcscf_object_init_ = true;
}

std::shared_ptr<SOMCSCF> CIWavefunction::mcscf_object(){
    if (!mcscf_object_init_){
        init_mcscf_object();
    }
    return somcscf_;
}


void CIWavefunction::print_vector(SharedCIVector vec, int root){
  int* mi_iac = init_int_array(Parameters_->nprint);
  int* mi_ibc = init_int_array(Parameters_->nprint);
  int* mi_iaidx = init_int_array(Parameters_->nprint);
  int* mi_ibidx = init_int_array(Parameters_->nprint);
  double* mi_coeff = init_array(Parameters_->nprint);

  // Print largest CI coefs
  vec->read(root, 0);
  vec->max_abs_vals(Parameters_->nprint, mi_iac, mi_ibc, mi_iaidx, mi_ibidx,
                    mi_coeff, Parameters_->neg_only);
  print_vec(Parameters_->nprint, mi_iac, mi_ibc, mi_iaidx, mi_ibidx, mi_coeff);

  free(mi_iac);
  free(mi_ibc);
  free(mi_iaidx);
  free(mi_ibidx);
  free(mi_coeff);
}

void CIWavefunction::compute_state_transfer(SharedCIVector ref, int ref_vec,
                                                SharedMatrix prop, SharedCIVector ret){

    if (!CalcInfo_->sigma_initialized){
        sigma_init(*(ref.get()), *(ret.get()));
    }

    // Figure out phase
    if (!Parameters_->Ms0) {
        int phase = 1;
    } else {
        int phase = ((int)Parameters_->S % 2) ? -1 : 1;
    }

    ref->read(ref_vec, 0);
    ret->zero();

    double** propp = prop->pointer();

    /* loop over unique I subblocks */
    // for (int Iblock = 0; Iblock < ref->num_blocks_; Iblock++) {
    //     // if (Parameters_->cc && !cc_reqd_Iblocks[Iblock]) continue;
    //     int Iacode = ref->Ia_code_[Iblock];
    //     int Ibcode = ref->Ib_code_[Iblock];
    //     int Inas = ref->Ia_size_[Iblock];
    //     int Inbs = ref->Ib_size_[Iblock];
    //     if (Inas == 0 || Inbs == 0) continue;
    //     double** coefs = ref->blocks_[Iblock];

    //     // Alpha components
    //     for (int Ia_idx = 0; Ia_idx < Inas; Ia_idx++){
    //         for (struct stringwr* Ib = betlist_[Ibcode], int Ib_idx = 0; Ib_idx < Inbs; Ib_idx++, Ib++){

    //         /* loop over excitations E^b_{ij} from |B(J_b)> */
    //         int Ibcnt = Ib->cnt[Iacode];
    //         unsigned int* Ibridx = Ib->ridx[Iacode];
    //         signed char* Ibsgn = Ib->sgn[Iacode];
    //         int* Iboij = Ib->oij[Iacode];
    //         for (Ib_ex = 0; Ib_ex < Ibcnt; Ib_ex++) {
    //           int oij = *Iboij++;
    //           int Ib_idx = *Ibridx++;
    //           int Ib_sgn = (double)*Ibsgn++;
    //           double c = ceoffs[Ia_idx][Ib_idx];
    //           int i = oij / CalcInfo_->num_ci_orbs;
    //           int j = oij % CalcInfo_->num_ci_orbs;
    //           propp[i][j] += c * c * Ib_sgn;
    //         }
    //     }


    //     } /* end loop over J blocks */
    // } /* end loop over J blocks */


    // if (ref->Ms0_) {
    //     if ((int)Parameters_->ref->% 2)
    //         ref->symmetrize(-1.0, 0);
    //     else
    //         ret->symmetrize(1.0, 0);
    // }

    // ret->write(ivec, 0);
}
}} // End Psi and CIWavefunction spaces
