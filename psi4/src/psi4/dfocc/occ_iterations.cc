/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2023 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This file is part of Psi4.
 *
 * Psi4 is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * Psi4 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along
 * with Psi4; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

#include "defines.h"
#include "dfocc.h"

#include "psi4/libqt/qt.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/writer.h"
#include "psi4/libmints/writer_file_prefix.h"
#include "psi4/libpsi4util/process.h"
#include "psi4/liboptions/liboptions.h"

using namespace psi;

namespace psi {
namespace dfoccwave {

void DFOCC::occ_iterations() {
    outfile->Printf("\n");
    outfile->Printf(" ============================================================================== \n");
    if (wfn_type_ == "DF-OMP2" && do_cd == "FALSE")
        outfile->Printf(" ================ Performing DF-OMP2 iterations... ============================ \n");
    else if (wfn_type_ == "DF-OMP3" && do_cd == "FALSE")
        outfile->Printf(" ================ Performing DF-OMP3 iterations... ============================ \n");
    else if (wfn_type_ == "DF-OMP2.5" && do_cd == "FALSE")
        outfile->Printf(" ================ Performing DF-OMP2.5 iterations... ========================== \n");
    else if (wfn_type_ == "DF-OLCCD" && do_cd == "FALSE")
        outfile->Printf(" ================ Performing DF-OLCCD iterations... =========================== \n");
    else if (wfn_type_ == "DF-OREMP" && do_cd == "FALSE")
        outfile->Printf(" ================ Performing DF-OREMP iterations... =========================== \n");
    else if (wfn_type_ == "DF-OMP2" && do_cd == "TRUE")
        outfile->Printf(" ================ Performing CD-OMP2 iterations... ============================ \n");
    else if (wfn_type_ == "DF-OMP3" && do_cd == "TRUE")
        outfile->Printf(" ================ Performing CD-OMP3 iterations... ============================ \n");
    else if (wfn_type_ == "DF-OMP2.5" && do_cd == "TRUE")
        outfile->Printf(" ================ Performing CD-OMP2.5 iterations... ============================ \n");
    else if (wfn_type_ == "DF-OLCCD" && do_cd == "TRUE")
        outfile->Printf(" ================ Performing CD-OLCCD iterations... ============================ \n");
    else if (wfn_type_ == "DF-OREMP" && do_cd == "TRUE")
        outfile->Printf(" ================ Performing CD-OREMP iterations... ============================ \n");
    outfile->Printf(" ============================================================================== \n");
    if (wfn_type_ == "DF-OMP2" && do_cd == "FALSE")
        outfile->Printf("\t            Minimizing DF-MP2-L Functional \n");
    else if (wfn_type_ == "DF-OMP3" && do_cd == "FALSE")
        outfile->Printf("\t            Minimizing DF-MP3-L Functional \n");
    else if (wfn_type_ == "DF-OMP2.5" && do_cd == "FALSE")
        outfile->Printf("\t            Minimizing DF-MP2.5-L Functional \n");
    else if (wfn_type_ == "DF-OLCCD" && do_cd == "FALSE")
        outfile->Printf("\t            Minimizing DF-LCCD-L Functional \n");
    else if (wfn_type_ == "DF-OREMP" && do_cd == "FALSE")
        outfile->Printf("\t            Minimizing DF-REMP-L Functional \n");
    else if (wfn_type_ == "DF-OMP2" && do_cd == "TRUE")
        outfile->Printf("\t            Minimizing CD-MP2-L Functional \n");
    else if (wfn_type_ == "DF-OMP3" && do_cd == "TRUE")
        outfile->Printf("\t            Minimizing CD-MP3-L Functional \n");
    else if (wfn_type_ == "DF-OMP2.5" && do_cd == "TRUE")
        outfile->Printf("\t            Minimizing CD-MP2.5-L Functional \n");
    else if (wfn_type_ == "DF-OLCCD" && do_cd == "TRUE")
        outfile->Printf("\t            Minimizing CD-LCCD-L Functional \n");
    else if (wfn_type_ == "DF-OREMP" && do_cd == "TRUE")
        outfile->Printf("\t            Minimizing CD-REMP-L Functional \n");
    outfile->Printf("\t            ------------------------------ \n");
    outfile->Printf(" Iter       E_total           DE           RMS MO Grad      MAX MO Grad      RMS T2    \n");
    outfile->Printf(" ----    ---------------    ----------     -----------      -----------     ---------- \n");

    //==========================================================================================
    //========================= NR iterations ==================================================
    //==========================================================================================
    itr_occ = 0;
    mu_ls = 0;
    conver = 1;  // Assuming that the MOs will be optimized.
    mo_optimized = 0;
    itr_diis = 0;
    ErempL = Emp2;

    // If diis?
    // if (noccA + noccB != 1) {
//    if (do_diis_ == 1) {
//        nvar = num_vecs + 1;
//        vecsA = std::make_shared<Tensor2d>("Alpha MO DIIS Vectors", num_vecs, nidpA);
//        errvecsA = std::make_shared<Tensor2d>("Alpha MO DIIS Error Vectors", num_vecs, nidpA);
//
//        if (reference_ == "UNRESTRICTED") {
//            vecsB = std::make_shared<Tensor2d>("Beta MO DIIS Vectors", num_vecs, nidpB);
//            errvecsB = std::make_shared<Tensor2d>("Beta MO DIIS Error Vectors", num_vecs, nidpB);
//        }
//    }
    //}

    //fire up DIIS
    if (do_diis_ == 1) {
//        outfile->Printf("firing up coupled DIIS...\n");
      orbitalDIIS = std::make_shared<DIISManager>(cc_maxdiis_, "Orbital Optimized DIIS", DIISManager::RemovalPolicy::LargestError, DIISManager::StoragePolicy::OnDisk); //initialize DIIS manager
      auto kappa_barA_ = std::make_shared<Vector>("Kappa_barA",nidpA);
      if (reference_ == "RESTRICTED" ) {
        if (wfn_type_ == "DF-OMP2" || wfn_type_ == "DF-OLCCD" ||  wfn_type_ == "DF-OREMP") {
          auto T2 = std::make_shared<Matrix>("T2", naoccA * navirA, naoccA * navirA); // T2 buffer prototype
          orbitalDIIS->set_error_vector_size(kappa_barA_.get(), T2.get());
          orbitalDIIS->set_vector_size(kappa_barA_.get(), T2.get());
          T2.reset();
        }
        else if (wfn_type_ == "DF-OMP2.5" || wfn_type_ == "DF-OMP3"){
          outfile->Printf("\ttrying to launch coupled DIIS for DF-OMP2.5/DF-OMP3...\n");
          auto T2 = std::make_shared<Matrix>("T2", naoccA * navirA, naoccA * navirA);   // T2_1 and T2_2 share the same prototype
          orbitalDIIS->set_error_vector_size(kappa_barA_.get(), T2.get(), T2.get());
      outfile->Printf("\tsuccessfully set the error vectors...\n");
          orbitalDIIS->set_vector_size(kappa_barA_.get(), T2.get(), T2.get());
      outfile->Printf("\tsuccessfully set the guess vectors...\n");
          T2.reset();
        }
//        else if (wfn_type_ == "DF-OLCCD" ||  wfn_type_ == "DF-OREMP"){
//         // still to come; PROBABLY IDENTICAL TO OMP2 -> save code...
//        }
      } else if (reference_ == "UNRESTRICTED"){
          auto kappa_barB_ = std::make_shared<Vector>("Kappa_barB", nidpB);
          if (wfn_type_ == "DF-OMP2"){
            auto T2AA = std::make_shared<Matrix>("T2AA", naoccA * navirA, naoccA * navirA);
            auto T2BB = std::make_shared<Matrix>("T2BB", naoccB * navirB, naoccB * navirB);
            auto T2AB = std::make_shared<Matrix>("T2AB", naoccA * navirA, naoccB * navirB);
            auto RT2AA = std::make_shared<Matrix>("RT2AA", naoccA * naoccA, navirA * navirA); // the amplitudes and the residuums differ in their signature
            auto RT2BB = std::make_shared<Matrix>("RT2BB", naoccB * naoccB, navirB * navirB); // a possible solution would be to reorder the residuum
            auto RT2AB = std::make_shared<Matrix>("RT2AB", naoccA * naoccB, navirA * navirB); // but this would be computationally wasted time
            orbitalDIIS->set_error_vector_size(kappa_barA_.get(), kappa_barB_.get(), RT2AA.get(), RT2BB.get(), RT2AB.get());
            orbitalDIIS->set_vector_size(kappa_barA_.get(), kappa_barB_.get(), T2AA.get(), T2BB.get(), T2AB.get());
            T2AA.reset();
            T2BB.reset();
            T2AB.reset();
            RT2AA.reset();
            RT2BB.reset();
            RT2AB.reset();
          }
          else if (wfn_type_ == "DF-OMP2.5" || wfn_type_ == "DF-OMP3"){
            auto T2AA = std::make_shared<Matrix>("T2AA", naoccA * naoccA, navirA * navirA); //T2_1 and T2_2 and their respective residuals share the same prototype
            auto T2BB = std::make_shared<Matrix>("T2BB", naoccB * naoccB, navirB * navirB);
            auto T2AB = std::make_shared<Matrix>("T2AB", naoccA * naoccB, navirA * navirB);
            orbitalDIIS->set_error_vector_size(kappa_barA_.get(), kappa_barB_.get(), T2AA.get(), T2BB.get(), T2AB.get(), T2AA.get(), T2BB.get(), T2AB.get());
            orbitalDIIS->set_vector_size(kappa_barA_.get(), kappa_barB_.get(), T2AA.get(), T2BB.get(), T2AB.get(), T2AA.get(), T2BB.get(), T2AB.get());
            T2AA.reset();
            T2BB.reset();
            T2AB.reset();
          }
          else if (wfn_type_ == "DF-OLCCD" || wfn_type_ == "DF-OREMP"){
            auto T2AA = std::make_shared<Matrix>("T2AA", naoccA * naoccA, navirA * navirA);
            auto T2BB = std::make_shared<Matrix>("T2BB", naoccB * naoccB, navirB * navirB);
            auto T2AB = std::make_shared<Matrix>("T2AB", naoccA * naoccB, navirA * navirB);
            orbitalDIIS->set_error_vector_size(kappa_barA_.get(), kappa_barB_.get(), T2AA.get(), T2BB.get(), T2AB.get());
            orbitalDIIS->set_vector_size(kappa_barA_.get(), kappa_barB_.get(), T2AA.get(), T2BB.get(), T2AB.get());
            T2AA.reset();
            T2BB.reset();
            T2AB.reset();
          }
          kappa_barB_.reset();
      }
      kappa_barA_.reset();
      outfile->Printf("\tsuccessfully fired up coupled DIIS...\n");
    }



    //==========================================================================================
    //========================= Head of the Loop ===============================================
    //==========================================================================================
    double last_rms_wog;
    do {
        itr_occ++;
        last_rms_wog = rms_wog;

        //==========================================================================================
        //========================= New orbital step ===============================================
        //==========================================================================================
        timer_on("kappa orb rot");
        if (hess_type == "HF") {
            if (orb_resp_solver_ == "LINEQ")
                kappa_orb_resp();
            else if (orb_resp_solver_ == "PCG")
                kappa_orb_resp_pcg();
        } else
            kappa_diag_hess();
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
        }  // end if (do_cd == "FALSE")

        // CD
        else if (do_cd == "TRUE") {
            timer_on("CD Trans");
            trans_cd();
            timer_off("CD Trans");
        }  // end if (do_cd == "TRUE")

        // Fock
        fock();

        // reference energy
        ref_energy();

        //==========================================================================================
        //========================= New Amplitudes =================================================
        //==========================================================================================
        if (wfn_type_ == "DF-OMP2")
            t2_1st_gen();

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

        else if (wfn_type_ == "DF-OREMP") {
            timer_on("T2");
            Fint_zero();
            remp_t2_amps();
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

        else if (wfn_type_ == "DF-OREMP") {
            lccd_pdm_3index_intr();
            omp3_opdm();
            oremp_tpdm();
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

        // else if (wfn_type_ == "DF-OMP3"  || wfn_type_ == "CD-OMP3") {
        else {
            gfock_cc_vo();
            gfock_cc_ov();
            gfock_cc_oo();
            gfock_cc_vv();
        }

        //==========================================================================================
        //========================= CCL ============================================================
        //==========================================================================================
        if (wfn_type_ == "DF-OMP2")
            mp2l_energy();
        else if (wfn_type_ == "DF-OMP3" || wfn_type_ == "DF-OMP2.5")
            mp3l_energy();
        else if (wfn_type_ == "DF-OLCCD" || wfn_type_ == "DF-OREMP")
            lccdl_energy();

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
            nidp = nidpA;
            rms_wog = rms_wogA;
            biggest_mograd = biggest_mogradA;
            rms_kappa = rms_kappaA;
            biggest_kappa = biggest_kappaA;
        }

        else if (reference_ == "UNRESTRICTED") {
            nidp = MAX0(nidpA, nidpB);
            rms_wog = MAX0(rms_wogA, rms_wogB);
            biggest_mograd = MAX0(biggest_mogradA, biggest_mogradB);
            rms_kappa = MAX0(rms_kappaA, rms_kappaB);
            biggest_kappa = MAX0(biggest_kappaA, biggest_kappaB);
            rms_t2 = MAX0(rms_t2AA, rms_t2BB);
            rms_t2 = MAX0(rms_t2, rms_t2AB);
        }

        if (wfn_type_ == "DF-OMP2") {
            outfile->Printf(" %3d     %12.10f  %12.2e   %12.2e     %12.2e    %12.2e \n", itr_occ, Emp2L, DE, rms_wog,
                            biggest_mograd, rms_t2);
        } else if (wfn_type_ == "DF-OMP3" || wfn_type_ == "DF-OMP2.5") {
            outfile->Printf(" %3d     %12.10f  %12.2e   %12.2e     %12.2e    %12.2e \n", itr_occ, Emp3L, DE, rms_wog,
                            biggest_mograd, rms_t2);
        } else if (wfn_type_ == "DF-OLCCD" || wfn_type_ == "DF-OREMP") {
            outfile->Printf(" %3d     %12.10f  %12.2e   %12.2e     %12.2e    %12.2e \n", itr_occ, ElccdL, DE, rms_wog,
                            biggest_mograd, rms_t2);
        }

        //==========================================================================================
        //========================= Convergence? ===================================================
        //==========================================================================================
        if (itr_occ >= mo_maxiter) {
            conver = 0;  // means MOs are NOT optimized
            break;
        }

        if (rms_wog >= DIVERGE) {
            throw PSIEXCEPTION("DF-OCC iterations are diverging");
        }

    } while (rms_wog >= tol_grad || biggest_mograd >= mograd_max || std::fabs(DE) >= tol_Eod || std::fabs(last_rms_wog - rms_wog) >= 10 * tol_Eod);

    if (conver == 1) {
        mo_optimized = 1;
        outfile->Printf("\n");
        outfile->Printf(" ============================================================================== \n");
        if (wfn_type_ == "DF-OMP2")
            outfile->Printf(" ======================== DF-OMP2 ITERATIONS ARE CONVERGED ==================== \n");
        else if (wfn_type_ == "DF-OMP3")
            outfile->Printf(" ======================== DF-OMP3 ITERATIONS ARE CONVERGED ==================== \n");
        else if (wfn_type_ == "DF-OMP2.5")
            outfile->Printf(" ======================== DF-OMP2.5 ITERATIONS ARE CONVERGED ================== \n");
        else if (wfn_type_ == "DF-OLCCD")
            outfile->Printf(" ======================== DF-OLCCD ITERATIONS ARE CONVERGED =================== \n");
        else if (wfn_type_ == "DF-OREMP")
            outfile->Printf(" ======================== DF-OREMP ITERATIONS ARE CONVERGED =================== \n");
        outfile->Printf(" ============================================================================== \n");

    }

    else if (conver == 0) {
        if (wfn_type_ == "DF-OMP2")
            outfile->Printf("\n ======================== DF-OMP2 IS NOT CONVERGED IN %2d ITERATIONS ========== \n",
                            mo_maxiter);
        else if (wfn_type_ == "DF-OMP3")
            outfile->Printf("\n ======================== DF-OMP3 IS NOT CONVERGED IN %2d ITERATIONS ========== \n",
                            mo_maxiter);
        else if (wfn_type_ == "DF-OMP2.5")
            outfile->Printf("\n ======================== DF-OMP2.5 IS NOT CONVERGED IN %2d ITERATIONS ======== \n",
                            mo_maxiter);
        else if (wfn_type_ == "DF-OLCCD")
            outfile->Printf("\n ======================== DF-OLCCD IS NOT CONVERGED IN %2d ITERATIONS ========= \n",
                            mo_maxiter);
        else if (wfn_type_ == "DF-OREMP")
            outfile->Printf("\n ======================== DF-OREMP IS NOT CONVERGED IN %2d ITERATIONS ========= \n",
                            mo_maxiter);

        throw PSIEXCEPTION("DF-OCC iterations did not converge");
    }

    if (do_diis_ == 1) {
      orbitalDIIS->delete_diis_file();
    }
}  // end occ_iterations

//=========================
// SAVE MOs to wfn
//=========================
void DFOCC::save_mo_to_wfn() {

    name_=wfn_type_.c_str();
    module_="dfocc";

    // Save MOs to wfn_; We cannot semicanonicalize them, as we'd need to do the same to all MO-basis quantities
    if (reference_ == "RESTRICTED") {
        SharedMatrix Ca = std::make_shared<Matrix>("Alpha MO Coefficients", nso_, nmo_);
        CmoA->to_shared_matrix(Ca);
        Ca_->copy(Ca);

        if (options_.get_str("MOLDEN_WRITE") == "TRUE") {
            // Diagonalize OPDM to obtain NOs
            SharedMatrix aevecs(new Matrix("Eigenvectors (Alpha)", nmo_, nmo_));
            SharedVector aevals(new Vector("Eigenvalues (Alpha)", nmo_));

            // Diagonaliz OPDM
            auto a_opdm = std::make_shared<Matrix>("Alpha OPDM", nmo_, nmo_);
            G1->to_shared_matrix(a_opdm);
            // scale by 1/2 because MoldenWrite expect only alpha part
            a_opdm->scale(0.5);
            a_opdm->diagonalize(aevecs, aevals, descending);

            // Form transformation matrix from AO to NO
            SharedMatrix aAONO(new Matrix("NOs (Alpha)", nso_, nmo_));
            aAONO->gemm(false, false, 1.0, Ca, aevecs, 0.0);

            // Write to MOLDEN file
            std::shared_ptr<MoldenWriter> molden(new MoldenWriter(shared_from_this()));
            std::string filename = get_writer_file_prefix(molecule_->name()) + "_dfocc.molden";

            // For now use zeros instead of energies, and NO occupation numbers as occupation numbers
            SharedVector dummy_a(new Vector("Dummy Vector Alpha", nmo_));
            for (int i = 0; i < naoccA; ++i) eps_orbA->set(i + nfrzc, eigooA->get(i));
            for (int a = 0; a < navirA; ++a) eps_orbA->set(a + noccA, eigvvA->get(a));
            eps_orbA->to_shared_vector(dummy_a);

            // write
            molden->write(filename, aAONO, aAONO, dummy_a, dummy_a, aevals, aevals, true);

            // free
            aAONO.reset();
            a_opdm.reset();
        }

        Ca.reset();
    }

    else if (reference_ == "UNRESTRICTED") {
        SharedMatrix Ca = std::make_shared<Matrix>("Alpha MO Coefficients", nso_, nmo_);
        SharedMatrix Cb = std::make_shared<Matrix>("Beta MO Coefficients", nso_, nmo_);
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
            SharedMatrix a_opdm = std::make_shared<Matrix>("Alpha OPDM", nmo_, nmo_);
            SharedMatrix b_opdm = std::make_shared<Matrix>("Alpha OPDM", nmo_, nmo_);
            G1A->to_shared_matrix(a_opdm);
            G1B->to_shared_matrix(b_opdm);
            a_opdm->diagonalize(aevecs, aevals, descending);
            b_opdm->diagonalize(bevecs, bevals, descending);

            // Form transformation matrix from AO to NO
            SharedMatrix aAONO(new Matrix("NOs (Alpha)", nso_, nmo_));
            SharedMatrix bAONO(new Matrix("NOs (Beta)", nso_, nmo_));
            aAONO->gemm(false, false, 1.0, Ca, aevecs, 0.0);
            bAONO->gemm(false, false, 1.0, Cb, bevecs, 0.0);

            // Write to MOLDEN file
            std::shared_ptr<MoldenWriter> molden(new MoldenWriter(shared_from_this()));
            std::string filename = get_writer_file_prefix(molecule_->name()) + "_dfocc.molden";

            // For now use zeros instead of energies, and NO occupation numbers as occupation numbers
            SharedVector dummy_a(new Vector("Dummy Vector Alpha", nmo_));
            SharedVector dummy_b(new Vector("Dummy Vector Beta", nmo_));
            for (int i = 0; i < naoccA; ++i) eps_orbA->set(i + nfrzc, eigooA->get(i));
            for (int a = 0; a < navirA; ++a) eps_orbA->set(a + noccA, eigvvA->get(a));
            for (int i = 0; i < naoccB; ++i) eps_orbB->set(i + nfrzc, eigooB->get(i));
            for (int a = 0; a < navirB; ++a) eps_orbB->set(a + noccB, eigvvB->get(a));
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

}  // end save_mo_to_wfn

// do a coupled DIIS extrapolation of the amplitudes and orbital rotations
// as is done in the canonical case
void DFOCC::oo_diis() {
  SharedTensor2d  t1, t2, t1aa, t1bb, t1ab, t2aa, t2bb, t2ab;
  SharedTensor2d  rt1, rt2, rt1aa, rt1bb, rt1ab, rt2aa, rt2bb, rt2ab;
  SharedTensor2d  U;

//  outfile->Printf("entering oo_diis...\n");
  if (!do_diis_) return;

    psio_->open(PSIF_DFOCC_AMPS, PSIO_OPEN_OLD);
    // current orbital rotations and orbital gradients
    auto wogA_vec = std::make_shared<Vector>("wogA",nidpA);
    wogA->to_shared_vector(wogA_vec);
    auto kappa_barA_vec = std::make_shared<Vector>("kappabarA",nidpA);
    kappa_barA->to_shared_vector(kappa_barA_vec);
    if (reference_ == "RESTRICTED"){
      if (wfn_type_ == "DF-OMP2")   {
        t1 = std::make_shared<Tensor2d>("T2_1 (ia|jb)", naoccA, navirA, naoccA, navirA);
        t1->read_symm(psio_, PSIF_DFOCC_AMPS);
        auto T1 = std::make_shared<Matrix>("T1",naoccA * navirA, naoccA * navirA);
        t1->to_matrix(T1);
        rt1 = std::make_shared<Tensor2d>("RT2_1 (ia|jb)", naoccA, navirA, naoccA, navirA);
        rt1->read_symm(psio_, PSIF_DFOCC_AMPS);
        auto RT1 = std::make_shared<Matrix>("RT1",naoccA * navirA, naoccA * navirA);
        rt1->to_matrix(RT1);
        orbitalDIIS->add_entry(wogA_vec.get(), RT1.get(), kappa_barA_vec.get(), T1.get());
        if (orbitalDIIS->subspace_size() >= cc_mindiis_) {
          orbitalDIIS->extrapolate(kappa_barA_vec.get(), T1.get());
          t1->set2(T1);
          t1->write_symm(psio_, PSIF_DFOCC_AMPS);
          t1.reset();
          for (int i=0; i<nidpA; i++ ){kappa_barA->set(i,kappa_barA_vec->get(i));}
        }
      }
      else if (wfn_type_ == "OMP2.5" || wfn_type_ == "OMP3") {
        t1 = std::make_shared<Tensor2d>("T2_1 (ia|jb)", naoccA, navirA, naoccA, navirA);
        t1->read_symm(psio_, PSIF_DFOCC_AMPS);
        auto T1 = std::make_shared<Matrix>("T1",naoccA * navirA, naoccA * navirA);
        t1->to_matrix(T1);
        t2 = std::make_shared<Tensor2d>("T2_2 (IA|JB)", naoccA, navirA, naoccA, navirA);
        t2->read_symm(psio_, PSIF_DFOCC_AMPS);
        auto T2 = std::make_shared<Matrix>("T2",naoccA * navirA, naoccA * navirA);
        t2->to_matrix(T2);
        rt1 = std::make_shared<Tensor2d>("RT2 (IA|JB)", naoccA, navirA, naoccA, navirA);
        rt1->read_symm(psio_, PSIF_DFOCC_AMPS);
        auto RT1 = std::make_shared<Matrix>("RT1",naoccA * navirA, naoccA * navirA);
        rt1->to_matrix(RT1);
        rt2 = std::make_shared<Tensor2d>("RT2_2 (IA|JB)", naoccA, navirA, naoccA, navirA);
        rt2->read_symm(psio_, PSIF_DFOCC_AMPS);
        auto RT2 = std::make_shared<Matrix>("RT2",naoccA * navirA, naoccA * navirA);
        rt2->to_matrix(RT2);
        orbitalDIIS->add_entry(wogA_vec.get(), RT1.get(), RT2.get(), kappa_barA_vec.get(), T1.get(), T2.get());
        if (orbitalDIIS->subspace_size() >= cc_mindiis_) {
          orbitalDIIS->extrapolate(kappa_barA_vec.get(), T1.get(), T2.get());
          t1->set2(T1);
          t1->write_symm(psio_, PSIF_DFOCC_AMPS);
          t2->set2(T2);
          t2->write_symm(psio_, PSIF_DFOCC_AMPS);
          for (int i=0; i<nidpA; i++ ){kappa_barA->set(i,kappa_barA_vec->get(i));}
        }
      }
      else if (wfn_type_ == "DF-OLCCD" || wfn_type_ == "DF-OREMP") {
        t1 = std::make_shared<Tensor2d>("T2 (IA|JB)", naoccA, navirA, naoccA, navirA);
        t1->read_symm(psio_, PSIF_DFOCC_AMPS);
        auto T1 = std::make_shared<Matrix>("T1",naoccA * navirA, naoccA * navirA);
        t1->to_matrix(T1);
        rt1 = std::make_shared<Tensor2d>("RT2 (IA|JB)", naoccA, navirA, naoccA, navirA);
        rt1->read_symm(psio_, PSIF_DFOCC_AMPS);
        auto RT1 = std::make_shared<Matrix>("RT1",naoccA * navirA, naoccA * navirA);
        rt1->to_matrix(RT1);
        orbitalDIIS->add_entry(wogA_vec.get(), RT1.get(), kappa_barA_vec.get(), T1.get());
        if (orbitalDIIS->subspace_size() >= cc_mindiis_) {
          orbitalDIIS->extrapolate(kappa_barA_vec.get(), T1.get());
          t1->set2(T1);
          t1->write_symm(psio_, PSIF_DFOCC_AMPS);
          for (int i=0; i<nidpA; i++ ){kappa_barA->set(i,kappa_barA_vec->get(i));}
        }
      }
    } else if (reference_ == "UNRESTRICTED") {
//      outfile->Printf("oo_diis:: entering unrestricted branch...\n");
        auto wogB_vec = std::make_shared<Vector>("wogB",nidpB);
        wogB->to_shared_vector(wogB_vec);
        auto kappa_barB_vec = std::make_shared<Vector>("kappabarB",nidpB);
        kappa_barB->to_shared_vector(kappa_barB_vec);

      if (wfn_type_ == "DF-OMP2")   {
         // AMPLITUDES
         // alpha-alpha
        t1aa = std::make_shared<Tensor2d>("T2_1 (IA|JB)", naoccA, navirA, naoccA, navirA);
        t1aa->read_symm(psio_, PSIF_DFOCC_AMPS);
//        outfile->Printf("read tha AA amplitudes\n");
        auto T1AA = std::make_shared<Matrix>("T1AA", naoccA * navirA, naoccA * navirA);
        t1aa->to_matrix(T1AA);
        // beta-beta
        t1bb = std::make_shared<Tensor2d>("T2_1 (ia|jb)", naoccB, navirB, naoccB, navirB);
        t1bb->read_symm(psio_, PSIF_DFOCC_AMPS);
//        outfile->Printf("read tha BB amplitudes\n");
        auto T1BB = std::make_shared<Matrix>("T1BB", naoccB * navirB, naoccB * navirB);
        t1bb->to_matrix(T1BB);
        // alpha-beta
        t1ab = std::make_shared<Tensor2d>("T2_1 (IA|jb)", naoccA, navirA, naoccB, navirB);
        t1ab->read(psio_, PSIF_DFOCC_AMPS);
//        outfile->Printf("read tha AB amplitudes\n");
        auto T1AB = std::make_shared<Matrix>("T1AB", naoccA * navirA, naoccB * navirB);
        t1ab->to_matrix(T1AB);
        // RESIDUALS
        // alpha-alpha
        rt1aa = std::make_shared<Tensor2d>("RT2_1 <IJ|AB>", naoccA, naoccA, navirA, navirA);
        rt1aa->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
//        outfile->Printf("read tha AA residuum\n");
        auto RT1AA = std::make_shared<Matrix>("RT1AA", naoccA * naoccA, navirA * navirA);
        rt1aa->to_matrix(RT1AA);
        // beta-beta
        rt1bb = std::make_shared<Tensor2d>("RT2_1 <ij|ab>", naoccB, naoccB, navirB, navirB);
        rt1bb->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
//        outfile->Printf("read tha BB residuum\n");
        auto RT1BB = std::make_shared<Matrix>("RT1BB", naoccB * naoccB, navirB * navirB);
        rt1bb->to_matrix(RT1BB);
        // alpha-beta
        rt1ab = std::make_shared<Tensor2d>("RT2_1 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
        rt1ab->read(psio_, PSIF_DFOCC_AMPS);
//        outfile->Printf("read tha AB residuum\n");
        auto RT1AB = std::make_shared<Matrix>("RT1AB", naoccA * naoccB, navirA * navirB);
        rt1ab->to_matrix(RT1AB);
//        outfile->Printf("successfully transferred the OMP2 amplitudes and residuals...\n");

//        outfile->Printf("successfully transferred the orbital gradient and rotation...\n");
        // all complete. now call add_entry, first pass the error vector componenets, then the guess vector components
        // in the same order as defined by set_error_vector_size and set_vector_size
        orbitalDIIS->add_entry(wogA_vec.get(), wogB_vec.get(), RT1AA.get(), RT1BB.get(), RT1AB.get(),
                               kappa_barA_vec.get(), kappa_barB_vec.get(), T1AA.get(), T1BB.get(), T1AB.get());
//        outfile->Printf("successfully added vectors to DIIS subspace\n");
        if (orbitalDIIS->subspace_size() >= cc_mindiis_) {
//          outfile->Printf("calling extrapolate... ");
          orbitalDIIS->extrapolate(kappa_barA_vec.get(), kappa_barB_vec.get(), T1AA.get(), T1BB.get(), T1AB.get());
//          outfile->Printf("success.\n ");
          t1aa->set2(T1AA);
          t1bb->set2(T1BB);
          t1ab->set2(T1AB);
          t1aa->write_symm(psio_, PSIF_DFOCC_AMPS);
//          t1aa.reset();
          t1bb->write_symm(psio_, PSIF_DFOCC_AMPS);
//          t1bb.reset();
          t1ab->write(psio_, PSIF_DFOCC_AMPS);
//          t1ab.reset();
          for (int i=0; i<nidpA; i++ ){kappa_barA->set(i,kappa_barA_vec->get(i));} // please don't tell me that this is the only way to get the rotations back to kappa_barA...
          for (int i=0; i<nidpB; i++ ){kappa_barB->set(i,kappa_barB_vec->get(i));} // if someone can do better: go for it
//          kappa_barA->set(&kappa_barA_vec);
//          kappa_barB->set(&kappa_barB_vec);

        }
          T1AA.reset();
          T1BB.reset();
          T1AB.reset();
          RT1AA.reset();
          RT1BB.reset();
          RT1AB.reset();
          kappa_barA_vec.reset();
          kappa_barB_vec.reset();
      } else if (wfn_type_ == "DF-OMP2.5" || wfn_type_ == "DF-OMP3") {
        // 1st order alpha-alpha amplitudes
          t1aa = std::make_shared<Tensor2d>("T2_1 <IJ|AB>", naoccA, naoccA, navirA, navirA);
          t1aa->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
          auto T1AA = std::make_shared<Matrix>("T1AA", naoccA * naoccA, navirA * navirA);
          t1aa->to_matrix(T1AA);
        // 2nd order alpha-alpha amplitudes
          t2aa = std::make_shared<Tensor2d>("T2_2 <IJ|AB>", naoccA, naoccA, navirA, navirA);
          t2aa->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
          auto T2AA = std::make_shared<Matrix>("T2AA", naoccA * naoccA, navirA * navirA);
          t2aa->to_matrix(T2AA);
        // 1st order beta-beta amplitudes
          t1bb = std::make_shared<Tensor2d>("T2_1 <ij|ab>", naoccB, naoccB, navirB, navirB);
          t1bb->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
          auto T1BB = std::make_shared<Matrix>("T1BB", naoccB * naoccB, navirB * navirB);
          t1bb->to_matrix(T1BB);
        // 2nd order beta amplitudes
          t2bb = std::make_shared<Tensor2d>("T2_2 <ij|ab>", naoccB, naoccB, navirB, navirB);
          t2bb->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
          auto T2BB = std::make_shared<Matrix>("T2BB", naoccB * naoccB, navirB * navirB);
          t2bb->to_matrix(T2BB);
        // 1st order alpha-beta amplitudes
          t1ab = std::make_shared<Tensor2d>("T2_1 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
          t1ab->read(psio_, PSIF_DFOCC_AMPS);
          auto T1AB = std::make_shared<Matrix>("T1AB", naoccA * naoccB, navirA * navirB);
          t1ab->to_matrix(T1AB);
        // 2nd order alpha-beta amplitudes
          t2ab = std::make_shared<Tensor2d>("T2_1 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
          t2ab->read(psio_, PSIF_DFOCC_AMPS);
          auto T2AB = std::make_shared<Matrix>("T2AB", naoccA * naoccB, navirA * navirB);
          t2ab->to_matrix(T2AB);
        // RESIDUALS
        // 1st order alpha-alpha residuals
          rt1aa = std::make_shared<Tensor2d>("RT2_1 <IJ|AB>", naoccA, naoccA, navirA, navirA);
          rt1aa->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
          auto RT1AA = std::make_shared<Matrix>("RT1AA", naoccA * naoccA, navirA * navirA);
          rt1aa->to_matrix(RT1AA);
        // 2nd order alpha-alpha residuals
          rt2aa = std::make_shared<Tensor2d>("RT2_2 <IJ|AB>", naoccA, naoccA, navirA, navirA);
          rt2aa->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
          auto RT2AA = std::make_shared<Matrix>("RT2AA", naoccA * naoccA, navirA * navirA);
          rt2aa->to_matrix(RT2AA);
        // 1st order beta-beta residuals
          rt1bb = std::make_shared<Tensor2d>("RT2_1 <ij|ab>", naoccB, naoccB, navirB, navirB);
          rt1bb->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
          auto RT1BB = std::make_shared<Matrix>("RT1BB", naoccB * naoccB, navirB * navirB);
          rt1bb->to_matrix(RT1BB);
        // 2nd order beta-beta residuals
          rt2bb = std::make_shared<Tensor2d>("RT2_2 <ij|ab>", naoccB, naoccB, navirB, navirB);
          rt2bb->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
          auto RT2BB = std::make_shared<Matrix>("RT2BB", naoccB * naoccB, navirB * navirB);
          rt2bb->to_matrix(RT2BB);
        // 1st order alpha-beta residuals
          rt1ab = std::make_shared<Tensor2d>("RT2_1 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
          rt1ab->read(psio_, PSIF_DFOCC_AMPS);
          auto RT1AB = std::make_shared<Matrix>("RT1AB", naoccA * naoccB, navirA * navirB);
          rt1ab->to_matrix(RT1AB);
       // 2nd order alpha-beta residuals
          rt2ab = std::make_shared<Tensor2d>("RT2_2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
          rt2ab->read(psio_, PSIF_DFOCC_AMPS);
          auto RT2AB = std::make_shared<Matrix>("RT2AB", naoccA * naoccB, navirA * navirB);
          rt2ab->to_matrix(RT2AB);
          orbitalDIIS->add_entry(wogA_vec.get(), wogB_vec.get(), RT1AA.get(), RT1BB.get(), RT1AB.get(),
                               RT2AA.get(), RT2BB.get(), RT2AB.get(),
                               kappa_barA_vec.get(), kappa_barB_vec.get(), T1AA.get(), T1BB.get(), T1AB.get(),
                               T2AA.get(), T2BB.get(), T2AB.get());
          if (orbitalDIIS->subspace_size() >= cc_mindiis_) {
            orbitalDIIS->extrapolate(kappa_barA_vec.get(), kappa_barB_vec.get(), T1AA.get(), T1BB.get(), T1AB.get(), T2AA.get(), T2BB.get(), T2AB.get());
            t1aa->set2(T1AA);
            t1bb->set2(T1BB);
            t1ab->set2(T1AB);
            t1aa->write_anti_symm(psio_, PSIF_DFOCC_AMPS);
            t1bb->write_anti_symm(psio_, PSIF_DFOCC_AMPS);
            t1ab->write(psio_, PSIF_DFOCC_AMPS);
            t2aa->set2(T2AA);
            t2bb->set2(T2BB);
            t2ab->set2(T2AB);
            t2aa->write_anti_symm(psio_, PSIF_DFOCC_AMPS);
            t2bb->write_anti_symm(psio_, PSIF_DFOCC_AMPS);
            t2ab->write(psio_, PSIF_DFOCC_AMPS);
            for (int i=0; i<nidpA; i++ ){kappa_barA->set(i,kappa_barA_vec->get(i));} // please don't tell me that this is the only way to get the rotations back to kappa_barA...
            for (int i=0; i<nidpB; i++ ){kappa_barB->set(i,kappa_barB_vec->get(i));} // if someone can do better: go for it

        }

      } else if (wfn_type_ == "DF-OLCCD" || wfn_type_ == "DF-OREMP") {
//          outfile->Printf("collecting data for DF-OLCCD DIIS\n");
          t1aa = std::make_shared<Tensor2d>("T2 <IJ|AB>", naoccA, naoccA, navirA, navirA);
          t1aa->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
          auto T1AA = std::make_shared<Matrix>("T1AA", naoccA * naoccA, navirA * navirA);
          t1aa->to_matrix(T1AA);
          t1bb = std::make_shared<Tensor2d>("T2 <ij|ab>", naoccB, naoccB, navirB, navirB);
          t1bb->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
          auto T1BB = std::make_shared<Matrix>("T1BB", naoccB * naoccB, navirB * navirB);
          t1bb->to_matrix(T1BB);
          t1ab = std::make_shared<Tensor2d>("T2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
          t1ab->read(psio_, PSIF_DFOCC_AMPS);
          auto T1AB = std::make_shared<Matrix>("T1AB", naoccA * naoccB, navirA * navirB);
          t1ab->to_matrix(T1AB);
          rt1aa = std::make_shared<Tensor2d>("RT2 <IJ|AB>", naoccA, naoccA, navirA, navirA);
          rt1aa->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
          auto RT1AA = std::make_shared<Matrix>("RT1AA", naoccA * naoccA, navirA * navirA);
          rt1aa->to_matrix(RT1AA);
          rt1bb = std::make_shared<Tensor2d>("RT2 <ij|ab>", naoccB, naoccB, navirB, navirB);
          rt1bb->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
          auto RT1BB = std::make_shared<Matrix>("RT1BB", naoccB * naoccB, navirB * navirB);
          rt1bb->to_matrix(RT1BB);
          rt1ab = std::make_shared<Tensor2d>("RT2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
          rt1ab->read(psio_, PSIF_DFOCC_AMPS);
          auto RT1AB = std::make_shared<Matrix>("RT1AB", naoccA * naoccB, navirA * navirB);
          rt1ab->to_matrix(RT1AB);
          orbitalDIIS->add_entry(wogA_vec.get(), wogB_vec.get(), RT1AA.get(), RT1BB.get(), RT1AB.get(),
                               kappa_barA_vec.get(), kappa_barB_vec.get(), T1AA.get(), T1BB.get(), T1AB.get());
          if (orbitalDIIS->subspace_size() >= cc_mindiis_) {
            orbitalDIIS->extrapolate(kappa_barA_vec.get(), kappa_barB_vec.get(), T1AA.get(), T1BB.get(), T1AB.get());
            t1aa->set2(T1AA);
            t1bb->set2(T1BB);
            t1ab->set2(T1AB);
            t1aa->write_anti_symm(psio_, PSIF_DFOCC_AMPS);
            t1bb->write_anti_symm(psio_, PSIF_DFOCC_AMPS);
            t1ab->write(psio_, PSIF_DFOCC_AMPS);
            // For some unknown reason, the OLCCD amplitudes are stored in two different formats in parallel
            // let's update both, otherwise the iterations will diverge.
            U = std::make_shared<Tensor2d>("T2 (IA|JB)", naoccA, navirA, naoccA, navirA);
            U->sort(1324, t1aa, 1.0, 0.0);
            U->write_symm(psio_, PSIF_DFOCC_AMPS);
            U.reset();
            U = std::make_shared<Tensor2d>("T2 (ia|ja)", naoccB, navirB, naoccB, navirB);
            U->sort(1324, t1bb, 1.0, 0.0);
            U->write_symm(psio_, PSIF_DFOCC_AMPS);
            U.reset();
            U = std::make_shared<Tensor2d>("T2 (IA|jb)", naoccA, navirA, naoccB, navirB);
            U->sort(1324, t1ab, 1.0, 0.0);
            U->write(psio_, PSIF_DFOCC_AMPS);
            U.reset();
            for (int i=0; i<nidpA; i++ ){kappa_barA->set(i,kappa_barA_vec->get(i));} // please don't tell me that this is the only way to get the rotations back to kappa_barA...
            for (int i=0; i<nidpB; i++ ){kappa_barB->set(i,kappa_barB_vec->get(i));} // if someone can do better: go for it


        }
      }
    }

    t1.reset();
    t2.reset();
    t1aa.reset();
    t1bb.reset();
    t1ab.reset();
    t2aa.reset();
    t2bb.reset();
    t2ab.reset();
    psio_->close(PSIF_DFOCC_AMPS,1);
} // end oo_diis


}  // namespace dfoccwave
}  // namespace psi
