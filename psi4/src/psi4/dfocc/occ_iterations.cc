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

#include "psi4/libqt/qt.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/writer.h"
#include "psi4/libmints/writer_file_prefix.h"
#include "defines.h"
#include "dfocc.h"

using namespace psi;
using namespace std;


namespace psi{ namespace dfoccwave{

void DFOCC::occ_iterations()
{

outfile->Printf("\n");
outfile->Printf(" ============================================================================== \n");
if (wfn_type_ == "DF-OMP2" && do_cd == "FALSE") outfile->Printf(" ================ Performing DF-OMP2 iterations... ============================ \n");
else if (wfn_type_ == "DF-OMP3" && do_cd == "FALSE") outfile->Printf(" ================ Performing DF-OMP3 iterations... ============================ \n");
else if (wfn_type_ == "DF-OMP2.5" && do_cd == "FALSE") outfile->Printf(" ================ Performing DF-OMP2.5 iterations... ========================== \n");
else if (wfn_type_ == "DF-OLCCD" && do_cd == "FALSE") outfile->Printf(" ================ Performing DF-OLCCD iterations... =========================== \n");
else if (wfn_type_ == "DF-OMP2" && do_cd == "TRUE") outfile->Printf(" ================ Performing CD-OMP2 iterations... ============================ \n");
else if (wfn_type_ == "DF-OMP3" && do_cd == "TRUE") outfile->Printf(" ================ Performing CD-OMP3 iterations... ============================ \n");
else if (wfn_type_ == "DF-OMP2.5" && do_cd == "TRUE") outfile->Printf(" ================ Performing CD-OMP2.5 iterations... ============================ \n");
else if (wfn_type_ == "DF-OLCCD" && do_cd == "TRUE") outfile->Printf(" ================ Performing CD-OLCCD iterations... ============================ \n");
outfile->Printf(" ============================================================================== \n");
if (wfn_type_ == "DF-OMP2" && do_cd == "FALSE") outfile->Printf( "\t            Minimizing DF-MP2-L Functional \n");
else if (wfn_type_ == "DF-OMP3" && do_cd == "FALSE") outfile->Printf( "\t            Minimizing DF-MP3-L Functional \n");
else if (wfn_type_ == "DF-OMP2.5" && do_cd == "FALSE") outfile->Printf( "\t            Minimizing DF-MP2.5-L Functional \n");
else if (wfn_type_ == "DF-OLCCD" && do_cd == "FALSE") outfile->Printf( "\t            Minimizing DF-LCCD-L Functional \n");
else if (wfn_type_ == "DF-OMP2" && do_cd == "TRUE") outfile->Printf( "\t            Minimizing CD-MP2-L Functional \n");
else if (wfn_type_ == "DF-OMP3" && do_cd == "TRUE") outfile->Printf( "\t            Minimizing CD-MP3-L Functional \n");
else if (wfn_type_ == "DF-OMP2.5" && do_cd == "TRUE") outfile->Printf( "\t            Minimizing CD-MP2.5-L Functional \n");
else if (wfn_type_ == "DF-OLCCD" && do_cd == "TRUE") outfile->Printf( "\t            Minimizing CD-LCCD-L Functional \n");
outfile->Printf( "\t            ------------------------------ \n");
outfile->Printf( " Iter       E_total           DE           RMS MO Grad      MAX MO Grad      RMS T2    \n");
outfile->Printf( " ----    ---------------    ----------     -----------      -----------     ---------- \n");



//==========================================================================================
//========================= NR iterations ==================================================
//==========================================================================================
      itr_occ = 0;
      mu_ls = 0;
      conver = 1; // Assuming that the MOs will be optimized.
      mo_optimized = 0;
      itr_diis = 0;

      // If diis?
      //if (noccA + noccB != 1) {
          if (do_diis_ == 1) {
              nvar = num_vecs +1;
              vecsA = SharedTensor2d(new Tensor2d("Alpha MO DIIS Vectors", num_vecs, nidpA));
              errvecsA = SharedTensor2d(new Tensor2d("Alpha MO DIIS Error Vectors", num_vecs, nidpA));

              if (reference_ == "UNRESTRICTED") {
                  vecsB = SharedTensor2d(new Tensor2d("Beta MO DIIS Vectors", num_vecs, nidpB));
                  errvecsB = SharedTensor2d(new Tensor2d("Beta MO DIIS Error Vectors", num_vecs, nidpB));
              }
          }
      //}

//==========================================================================================
//========================= Head of the Loop ===============================================
//==========================================================================================
do
{
       itr_occ++;

//==========================================================================================
//========================= New orbital step ===============================================
//==========================================================================================
        timer_on("kappa orb rot");
        if (hess_type == "HF") {
           if (orb_resp_solver_ == "LINEQ") kappa_orb_resp();
           else if (orb_resp_solver_ == "PCG") kappa_orb_resp_pcg();
        }
        else kappa_diag_hess();
        timer_off("kappa orb rot");

//==========================================================================================
//========================= update mo coefficients =========================================
//==========================================================================================
        timer_on("update_mo");
        update_mo();
        timer_off("update_mo");

//==========================================================================================
//========================= Trans TEI ======================================================
//==========================================================================================
    // DF
    if (do_cd == "FALSE") {
        timer_on("DF CC Integrals");
        trans_corr();
        timer_off("DF CC Integrals");

        timer_on("DF REF Integrals");
        trans_ref();
        timer_off("DF REF Integrals");
    }// end if (do_cd == "FALSE")

    // CD
    else if (do_cd == "TRUE") {
        timer_on("CD Trans");
        trans_cd();
        timer_off("CD Trans");
    }// end if (do_cd == "TRUE")

        // Fock
        fock();

        // reference energy
        ref_energy();

//==========================================================================================
//========================= New Amplitudes =================================================
//==========================================================================================
     if (wfn_type_ == "DF-OMP2") t2_1st_gen();

     else if (wfn_type_ == "DF-OMP3" || wfn_type_ == "DF-OMP2.5") {
	mp3_t2_1st_gen();
        timer_on("T2(2)");
	t2_2nd_gen();
        timer_off("T2(2)");
     }

     else if (wfn_type_ == "DF-OLCCD") {
        timer_on("T2");
	Fint_zero();
	lccd_t2_amps();
        timer_off("T2");
     }

//==========================================================================================
//========================= PDMs ===========================================================
//==========================================================================================
	if (wfn_type_ == "DF-OMP2") {
	    omp2_opdm();
	    omp2_tpdm();
            separable_tpdm();
        }

        else if (wfn_type_ == "DF-OMP3" || wfn_type_ == "DF-OMP2.5") {
	    mp3_pdm_3index_intr();
	    omp3_opdm();
	    omp3_tpdm();
	    sep_tpdm_cc();
        }

        else if (wfn_type_ == "DF-OLCCD") {
	    lccd_pdm_3index_intr();
	    omp3_opdm();
	    olccd_tpdm();
	    sep_tpdm_cc();
        }

//==========================================================================================
//========================= GFM ============================================================
//==========================================================================================
	if (wfn_type_ == "DF-OMP2") {
            gfock_vo();
            gfock_ov();
            gfock_oo();
            gfock_vv();
        }

	//else if (wfn_type_ == "DF-OMP3"  || wfn_type_ == "CD-OMP3") {
	else {
	    gfock_cc_vo();
	    gfock_cc_ov();
	    gfock_cc_oo();
	    gfock_cc_vv();
        }

//==========================================================================================
//========================= CCL ============================================================
//==========================================================================================
        if (wfn_type_ == "DF-OMP2") mp2l_energy();
        else if (wfn_type_ == "DF-OMP3" || wfn_type_ == "DF-OMP2.5") mp3l_energy();
        else if (wfn_type_ == "DF-OLCCD") lccdl_energy();


//==========================================================================================
//========================= MO Grad ========================================================
//==========================================================================================
        timer_on("MO Grad");
	mograd();
        timer_off("MO Grad");

//==========================================================================================
//========================= Print ==========================================================
//==========================================================================================
    if (reference_ == "RESTRICTED") {
	nidp=nidpA;
	rms_wog=rms_wogA;
	biggest_mograd=biggest_mogradA;
	rms_kappa=rms_kappaA;
	biggest_kappa=biggest_kappaA;
    }

    else if (reference_ == "UNRESTRICTED") {
	nidp=MAX0(nidpA,nidpB);
	rms_wog=MAX0(rms_wogA,rms_wogB);
	biggest_mograd=MAX0(biggest_mogradA,biggest_mogradB);
	rms_kappa=MAX0(rms_kappaA,rms_kappaB);
	biggest_kappa=MAX0(biggest_kappaA,biggest_kappaB);
	rms_t2=MAX0(rms_t2AA,rms_t2BB);
	rms_t2=MAX0(rms_t2,rms_t2AB);
    }

if(wfn_type_ == "DF-OMP2") {
	outfile->Printf(" %3d     %12.10f  %12.2e   %12.2e     %12.2e    %12.2e \n",itr_occ,Emp2L,DE,rms_wog,biggest_mograd,rms_t2);
}
else if(wfn_type_ == "DF-OMP3" || wfn_type_ == "DF-OMP2.5") {
	outfile->Printf(" %3d     %12.10f  %12.2e   %12.2e     %12.2e    %12.2e \n",itr_occ,Emp3L,DE,rms_wog,biggest_mograd,rms_t2);
}
else if(wfn_type_ == "DF-OLCCD") outfile->Printf(" %3d     %12.10f  %12.2e   %12.2e     %12.2e    %12.2e \n",itr_occ,ElccdL,DE,rms_wog,biggest_mograd,rms_t2);


//==========================================================================================
//========================= Convergence? ===================================================
//==========================================================================================
    if (itr_occ >= mo_maxiter) {
      conver = 0; // means MOs are NOT optimized
      break;
    }

    if (wfn_type_ != "DF-OLCCD") {
        if (rms_wog < tol_grad && biggest_mograd < mograd_max) break;
        if (fabs(DE) <= tol_Eod) break;
    }

    if (rms_wog >= DIVERGE) {
        throw PSIEXCEPTION("DF-OCC iterations are diverging");
    }

}
while(rms_wog >= tol_grad || biggest_mograd >= mograd_max || fabs(DE) >= tol_Eod);

if (conver == 1) {
mo_optimized = 1;
outfile->Printf("\n");
outfile->Printf(" ============================================================================== \n");
if (wfn_type_ == "DF-OMP2") outfile->Printf(" ======================== DF-OMP2 ITERATIONS ARE CONVERGED ==================== \n");
else if (wfn_type_ == "DF-OMP3") outfile->Printf(" ======================== DF-OMP3 ITERATIONS ARE CONVERGED ==================== \n");
else if (wfn_type_ == "DF-OMP2.5") outfile->Printf(" ======================== DF-OMP2.5 ITERATIONS ARE CONVERGED ================== \n");
else if (wfn_type_ == "DF-OLCCD") outfile->Printf(" ======================== DF-OLCCD ITERATIONS ARE CONVERGED =================== \n");
outfile->Printf(" ============================================================================== \n");

}

else if (conver == 0) {
  if (wfn_type_ == "DF-OMP2") outfile->Printf("\n ======================== DF-OMP2 IS NOT CONVERGED IN %2d ITERATIONS ========== \n", mo_maxiter);
  else if (wfn_type_ == "DF-OMP3") outfile->Printf("\n ======================== DF-OMP3 IS NOT CONVERGED IN %2d ITERATIONS ========== \n", mo_maxiter);
  else if (wfn_type_ == "DF-OMP2.5") outfile->Printf("\n ======================== DF-OMP2.5 IS NOT CONVERGED IN %2d ITERATIONS ======== \n", mo_maxiter);
  else if (wfn_type_ == "DF-OLCCD") outfile->Printf("\n ======================== DF-OLCCD IS NOT CONVERGED IN %2d ITERATIONS ========= \n", mo_maxiter);

  throw PSIEXCEPTION("DF-OCC iterations did not converge");
}

}// end occ_iterations


//=========================
// SAVE MOs to wfn
//=========================
void DFOCC::save_mo_to_wfn()
{
    // make sure we have semicanonic MOs
    if (orbs_already_sc == 0) semi_canonic();

    // Save mos to wfn
    if (reference_ == "RESTRICTED") {
	SharedMatrix Ca = SharedMatrix(new Matrix("Alpha MO Coefficients", nso_, nmo_));
	CmoA->to_shared_matrix(Ca);
    Ca_->copy(Ca);

      if (options_.get_str("MOLDEN_WRITE") == "TRUE") {
	// Diagonalize OPDM to obtain NOs
	SharedMatrix aevecs(new Matrix("Eigenvectors (Alpha)", nmo_, nmo_));
	SharedVector aevals(new Vector("Eigenvalues (Alpha)", nmo_));

	// Diagonaliz OPDM
	SharedMatrix a_opdm = SharedMatrix(new Matrix("Alpha OPDM", nmo_, nmo_));
	G1->to_shared_matrix(a_opdm);
	// scale by 1/2 because MoldenWrite expect only alpha part
	a_opdm->scale(0.5);
        a_opdm->diagonalize(aevecs, aevals, descending);

	// Form transformation matrix from AO to NO
        SharedMatrix aAONO (new Matrix("NOs (Alpha)", nso_, nmo_));
	aAONO->gemm(false, false, 1.0, Ca, aevecs, 0.0);

	// Write to MOLDEN file
	std::shared_ptr<MoldenWriter> molden(new MoldenWriter(shared_from_this()));
	std::string filename = get_writer_file_prefix(molecule_->name()) + "_dfocc.molden";

        // For now use zeros instead of energies, and DCFT NO occupation numbers as occupation numbers
	SharedVector dummy_a(new Vector("Dummy Vector Alpha", nmo_));
        for(int i = 0; i < naoccA; ++i) eps_orbA->set(i+nfrzc, eigooA->get(i));
        for(int a = 0; a < navirA; ++a) eps_orbA->set(a+noccA, eigvvA->get(a));
	eps_orbA->to_shared_vector(dummy_a);

	// write
	molden->write(filename, aAONO, aAONO, dummy_a, dummy_a, aevals, aevals, true);

	//free
	aAONO.reset();
	a_opdm.reset();
      }

	Ca.reset();
    }

    else if (reference_ == "UNRESTRICTED") {
	SharedMatrix Ca = SharedMatrix(new Matrix("Alpha MO Coefficients", nso_, nmo_));
	SharedMatrix Cb = SharedMatrix(new Matrix("Beta MO Coefficients", nso_, nmo_));
	CmoA->to_shared_matrix(Ca);
	CmoB->to_shared_matrix(Cb);

    Ca_->copy(Ca);
    Cb_->copy(Cb);

      if (options_.get_str("MOLDEN_WRITE") == "TRUE") {
	// Diagonalize OPDM to obtain NOs
	SharedMatrix aevecs(new Matrix("Eigenvectors (Alpha)", nmo_, nmo_));
	SharedMatrix bevecs(new Matrix("Eigenvectors (Beta)", nmo_, nmo_));
	SharedVector aevals(new Vector("Eigenvalues (Alpha)", nmo_));
	SharedVector bevals(new Vector("Eigenvalues (Beta)", nmo_));

	// Diagonaliz OPDM
	SharedMatrix a_opdm = SharedMatrix(new Matrix("Alpha OPDM", nmo_, nmo_));
	SharedMatrix b_opdm = SharedMatrix(new Matrix("Alpha OPDM", nmo_, nmo_));
	G1A->to_shared_matrix(a_opdm);
	G1B->to_shared_matrix(b_opdm);
        a_opdm->diagonalize(aevecs, aevals, descending);
        b_opdm->diagonalize(bevecs, bevals, descending);

	// Form transformation matrix from AO to NO
        SharedMatrix aAONO (new Matrix("NOs (Alpha)", nso_, nmo_));
        SharedMatrix bAONO (new Matrix("NOs (Beta)", nso_, nmo_));
	aAONO->gemm(false, false, 1.0, Ca, aevecs, 0.0);
	bAONO->gemm(false, false, 1.0, Cb, bevecs, 0.0);

	// Write to MOLDEN file
	std::shared_ptr<MoldenWriter> molden(new MoldenWriter(shared_from_this()));
	std::string filename = get_writer_file_prefix(molecule_->name()) + "_dfocc.molden";

        // For now use zeros instead of energies, and DCFT NO occupation numbers as occupation numbers
	SharedVector dummy_a(new Vector("Dummy Vector Alpha", nmo_));
	SharedVector dummy_b(new Vector("Dummy Vector Beta", nmo_));
        for(int i = 0; i < naoccA; ++i) eps_orbA->set(i+nfrzc, eigooA->get(i));
        for(int a = 0; a < navirA; ++a) eps_orbA->set(a+noccA, eigvvA->get(a));
        for(int i = 0; i < naoccB; ++i) eps_orbB->set(i+nfrzc, eigooB->get(i));
        for(int a = 0; a < navirB; ++a) eps_orbB->set(a+noccB, eigvvB->get(a));
	eps_orbA->to_shared_vector(dummy_a);
	eps_orbB->to_shared_vector(dummy_b);

	// write
	molden->write(filename, aAONO, bAONO, dummy_a, dummy_b, aevals, bevals, true);

	// free
	aAONO.reset();
	bAONO.reset();
	a_opdm.reset();
	b_opdm.reset();
      }

	Ca.reset();
	Cb.reset();
    }

} // end save_mo_to_wfn

}} // End Namespaces
