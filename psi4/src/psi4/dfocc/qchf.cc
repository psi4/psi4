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
#include "defines.h"
#include "dfocc.h"
#include "psi4/libciomr/libciomr.h"

using namespace psi;
using namespace std;


namespace psi{ namespace dfoccwave{

void DFOCC::qchf()
{
outfile->Printf("\n");
outfile->Printf(" ============================================================================== \n");
outfile->Printf(" ================ Performing QCHF iterations... =============================== \n");
outfile->Printf(" ============================================================================== \n");
outfile->Printf( "\t            QCHF \n");
outfile->Printf( "\t   ------------------------------ \n");
outfile->Printf( " Iter       E_total           DE           RMS MO Grad      MAX MO Grad  \n");
outfile->Printf( " ----    ---------------    ----------     -----------      -----------  \n");

//==========================================================================================
//========================= NR iterations ==================================================
//==========================================================================================
      itr_occ = 0;
      mu_ls = 0;
      int conver_hf = 1; // Assuming that the MOs will be optimized.
      double Eref_old = 0.0;

      // gwh
      gwh();

    // DF
    if (do_cd == "FALSE") {
        timer_on("DF REF Integrals");
        trans_ref();
        timer_off("DF REF Integrals");
        // Trans OEI
        timer_on("Trans OEI");
        trans_oei();
        timer_off("Trans OEI");
    }// end if (do_cd == "FALSE")

    // CD
    else if (do_cd == "TRUE") {
        timer_on("CD Trans");
        trans_cd();
        timer_off("CD Trans");
    }// end if (do_cd == "TRUE")

      // fock
      fock();
      //CmoA->print();

      // idp
      idp_hf();

     // MO Grad
     if (reference_ == "RESTRICTED") {
         WorbA->copy(FockA);
         WorbA->scale(2.0);
     }
     else if (reference_ == "UNRESTRICTED") {
         WorbA->copy(FockA);
         WorbB->copy(FockB);
         WorbA->scale(2.0);
         WorbB->scale(2.0);
     }

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
        kappa_qchf();
        timer_off("kappa orb rot");

//==========================================================================================
//========================= update mo coefficients =========================================
//==========================================================================================
        timer_on("update_mo");
        update_hfmo();
        timer_off("update_mo");

//==========================================================================================
//========================= Trans TEI ======================================================
//==========================================================================================
    // DF
    if (do_cd == "FALSE") {
        timer_on("DF REF Integrals");
        trans_ref();
        timer_off("DF REF Integrals");
        // Trans OEI
        timer_on("Trans OEI");
        trans_oei();
        timer_off("Trans OEI");
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
        Escf = Eref;
        DE = Eref - Eref_old;
        Eref_old = Eref;

//==========================================================================================
//========================= MO Grad ========================================================
//==========================================================================================
        timer_on("MO Grad");
  if (reference_ == "RESTRICTED") {
      WorbA->copy(FockA);
      WorbA->scale(2.0);
      //rms_wogA = 2.0*FvoA->rms();

      // Form w vector
      // Alpha
      for(int x = 0; x < nidpA; x++) {
	int p = idprowA->get(x);
	int q = idpcolA->get(x);
	wogA->set(x, WorbA->get(p,q));
      }
      //wogA->print();

    // find biggest_mograd
    biggest_mogradA=0;
    for (int i=0; i<nidpA;i++){
      if (fabs(wogA->get(i)) > biggest_mogradA)  biggest_mogradA=fabs(wogA->get(i));
    }

    // rms
    rms_wogA=0;
    for (int i=0; i<nidpA;i++) rms_wogA += wogA->get(i) * wogA->get(i);
    rms_wogA = wogA->rms();
    rms_wog=rms_wogA;

  }// end RHF

  else if (reference_ == "UNRESTRICTED") {
      WorbA->copy(FockA);
      WorbB->copy(FockB);
      WorbA->scale(2.0);
      WorbB->scale(2.0);
      //rms_wogA = 2.0*FvoA->rms();
      //rms_wogB = 2.0*FvoB->rms();

      // Form w vector
      // Alpha
      for(int x = 0; x < nidpA; x++) {
	int p = idprowA->get(x);
	int q = idpcolA->get(x);
	wogA->set(x, WorbA->get(p,q));
      }
      //wogA->print();

      // Beta
      for(int x = 0; x < nidpB; x++) {
	int p = idprowB->get(x);
	int q = idpcolB->get(x);
	wogB->set(x, WorbB->get(p,q));
      }
      //wogB->print();

    // find biggest_mograd
    biggest_mogradA=0;
    for (int i=0; i<nidpA;i++){
      if (fabs(wogA->get(i)) > biggest_mogradA)  biggest_mogradA=fabs(wogA->get(i));
    }

    biggest_mogradB=0;
    for (int i=0; i<nidpB;i++){
      if (fabs(wogB->get(i)) > biggest_mogradB)  biggest_mogradB=fabs(wogB->get(i));
    }

    // rms
    rms_wogA=0;
    for (int i=0; i<nidpA;i++) rms_wogA += wogA->get(i) * wogA->get(i);
    rms_wogA = wogA->rms();

    rms_wogB=0;
    for (int i=0; i<nidpB;i++) rms_wogB += wogB->get(i) * wogB->get(i);
    rms_wogB = wogB->rms();
    rms_wog=MAX0(rms_wogA,rms_wogB);
  }// end UHF
        timer_off("MO Grad");

//==========================================================================================
//========================= Print ==========================================================
//==========================================================================================
    if (reference_ == "RESTRICTED") {
	nidp=nidpA;
	biggest_mograd=biggest_mogradA;
    }

    else if (reference_ == "UNRESTRICTED") {
	nidp=MAX0(nidpA,nidpB);
	biggest_mograd=MAX0(biggest_mogradA,biggest_mogradB);
    }

outfile->Printf(" %3d     %12.10f  %12.2e   %12.2e     %12.2e \n",itr_occ,Eref,DE,rms_wog,biggest_mograd);

//==========================================================================================
//========================= Convergence? ===================================================
//==========================================================================================
    if (itr_occ >= mo_maxiter) {
      conver = 0; // means MOs are NOT optimized
      break;
    }

    if (rms_wog < tol_grad && biggest_mograd < mograd_max) break;
    if (fabs(DE) <= tol_Eod) break;

    if (rms_wog >= DIVERGE) {
        throw PSIEXCEPTION("QCHF iterations are diverging");
    }

}
while(rms_wog >= tol_grad || biggest_mograd >= mograd_max);

if (conver_hf == 1) {
outfile->Printf("\n");
outfile->Printf(" ============================================================================== \n");
outfile->Printf(" ======================== QCHF ITERATIONS ARE CONVERGED ======================= \n");
outfile->Printf(" ============================================================================== \n");


    // canonicalize
    canonic();

    // DF
    if (do_cd == "FALSE") {
        timer_on("DF REF Integrals");
        trans_ref();
        timer_off("DF REF Integrals");
        // Trans OEI
        timer_on("Trans OEI");
        trans_oei();
        timer_off("Trans OEI");
    }// end if (do_cd == "FALSE")

    // CD
    else if (do_cd == "TRUE") {
        timer_on("CD Trans");
        trans_cd();
        timer_off("CD Trans");
    }// end if (do_cd == "TRUE")

    // Fock
    fock();

    // memfree
    wogA.reset();
    kappaA.reset();
    kappa_newA.reset();
    kappa_barA.reset();
    idprowA.reset();
    idpcolA.reset();
    if (reference_ == "UNRESTRICTED") {
        wogB.reset();
        kappaB.reset();
        kappa_newB.reset();
        kappa_barB.reset();
        idprowB.reset();
        idpcolB.reset();
    }

}// end conver-hf = 1

else if (conver_hf == 0) {
  outfile->Printf("\n\tQCHF IS NOT CONVERGED IN %2d ITERATIONS ========== \n", mo_maxiter);
  throw PSIEXCEPTION("QCHF iterations did not converge");
}

}// end qchf

//=======================================================
//          IDP (HF)
//=======================================================
void DFOCC::idp_hf()
{
     int dim, block;

if (reference_ == "RESTRICTED") {
    // Form IDPs
    nidpA=0;

    // All V-O
    nidpA += nvirA * noccA;

    if (nidpA > 0) {
      idp_returnA = 1;
      wogA = SharedTensor1d(new Tensor1d("Alpha MO grad vector", nidpA));
      kappaA = SharedTensor1d(new Tensor1d("Alpha orb rot params vector of current step", nidpA));
      kappa_newA = SharedTensor1d(new Tensor1d("Alpha new orb rot params vector of current step", nidpA));
      kappa_barA = SharedTensor1d(new Tensor1d("Alpha orb rot params vector with respect to scf MOs", nidpA));
      //wog_intA = SharedTensor1d(new Tensor1d("Alpha Interpolated MO grad vector", nidpA));
      idprowA = SharedTensor1i(new Tensor1i("Alpha IDP Row", nidpA));
      idpcolA = SharedTensor1i(new Tensor1i("Alpha IDP Col", nidpA));

      // set idpA
      dim=0;

      // V-O
      //if (nfrzc == 0 && nfrzv == 0) {
          for(int a = 0; a < nvirA; a++){
	      for(int i = 0; i < noccA; i++){
	          idprowA->set(dim, a + noccA);
	          idpcolA->set(dim, i);
	          dim++;
	      }
          }
      //}
    }// end if nidpA != 0

    else if (nidpA == 0) {
            outfile->Printf("\tThere is not any non-redundant orbital rotation pair! \n");
            tstop();
            exit(EXIT_SUCCESS);
    }

}// end if (reference_ == "RESTRICTED")

else if (reference_ == "UNRESTRICTED") {
    // Form IDPs
    nidpA=0;
    nidpB=0;

    // All V-O
    nidpA += nvirA * noccA;
    nidpB += nvirB * noccB;

    if (nidpA == 0 && nidpB == 0) {
        outfile->Printf("\tThere is not any non-redundant orbital rotation pair! \n");
        tstop();
        exit(EXIT_SUCCESS);
    }

     if (nidpA > 0) {
      idp_returnA = 1;
      wogA = SharedTensor1d(new Tensor1d("Alpha MO grad vector", nidpA));
      kappaA = SharedTensor1d(new Tensor1d("Alpha orb rot params vector of current step", nidpA));
      kappa_newA = SharedTensor1d(new Tensor1d("Alpha new orb rot params vector of current step", nidpA));
      kappa_barA = SharedTensor1d(new Tensor1d("Alpha orb rot params vector with respect to scf MOs", nidpA));
      //wog_intA = SharedTensor1d(new Tensor1d("Alpha Interpolated MO grad vector", nidpA));
      idprowA = SharedTensor1i(new Tensor1i("Alpha IDP Row", nidpA));
      idpcolA = SharedTensor1i(new Tensor1i("Alpha IDP Col", nidpA));

      // set idpA
      dim=0;

      // V-O
      //if (nfrzc == 0 && nfrzv == 0) {
          for(int a = 0; a < nvirA; a++){
	      for(int i = 0; i < noccA; i++){
	          idprowA->set(dim, a + noccA);
	          idpcolA->set(dim, i);
	          dim++;
	      }
          }
      //}
    }// end if nidpA != 0

     if (nidpB > 0) {
      idp_returnB = 1;
      wogB = SharedTensor1d(new Tensor1d("Beta MO grad vector", nidpB));
      kappaB = SharedTensor1d(new Tensor1d("Beta orb rot params vector of current step", nidpB));
      kappa_newB = SharedTensor1d(new Tensor1d("Beta new orb rot params vector of current step", nidpB));
      kappa_barB = SharedTensor1d(new Tensor1d("Beta orb rot params vector with respect to scf MOs", nidpB));
      //wog_intB = SharedTensor1d(new Tensor1d("Beta Interpolated MO grad vector", nidpB));
      idprowB = SharedTensor1i(new Tensor1i("Beta IDP Row", nidpB));
      idpcolB = SharedTensor1i(new Tensor1i("Beta IDP Col", nidpB));

      // set idpB
      dim=0;

      // V-O
      //if (nfrzc == 0 && nfrzv == 0) {
          for(int a = 0; a < nvirB; a++){
	      for(int i = 0; i < noccB; i++){
	          idprowB->set(dim, a + noccB);
	          idpcolB->set(dim, i);
	          dim++;
	      }
          }
      //}
    }// end if nidpB != 0

}// end if (reference_ == "UNRESTRICTED")
}// end of idp_hf

//=======================================================
//          GWH
//=======================================================
void DFOCC::gwh()
{
     // Memalloc
     SharedTensor2d Fso = SharedTensor2d(new Tensor2d("SO-basis Fock Matrix", nso_, nso_));
     SharedTensor2d Fsop = SharedTensor2d(new Tensor2d("SO-basis Fock' Matrix", nso_, nso_));
     SharedTensor2d Smhalf = SharedTensor2d(new Tensor2d("S^-1/2", nso_, nso_));
     SharedTensor2d Cmop = SharedTensor2d(new Tensor2d("C' matrix", nso_, nmo_));
     SharedTensor2d Uso = SharedTensor2d(new Tensor2d("SO-basis U", nso_, nso_));
     SharedTensor2d temp = SharedTensor2d(new Tensor2d("Temp", nso_, nso_));
     SharedTensor1d e_orb = std::shared_ptr<Tensor1d>(new Tensor1d("epsilon <n|n>", nso_));
     SharedTensor1d DiagS = std::shared_ptr<Tensor1d>(new Tensor1d("Diag S", nso_));

     // F_mn = 1/2 * S_mn (H_mm + H_nn)
     for (int mu = 0; mu < nso_; mu++){
          for (int nu = 0; nu < nso_; nu++){
               double value = 0.875 * Sso->get(mu, nu) * ( Hso->get(mu, mu) + Hso->get(nu, nu) );
               Fso->set(mu, nu, value);
          }
     }

     // Diagonalize
     Sso->diagonalize(Uso, DiagS, cutoff);

     // Form S^(-1/2)
     for (int p = 0; p < nso_; p++) {
          DiagS->set(p, 1/sqrt(DiagS->get(p)));
     }

     for (int p = 0; p < nso_; p++) {
          Smhalf->set(p, p, DiagS->get(p));
     }

     // Diagonalize Fock matrix
     temp->gemm(true, false, Smhalf, Fso, 1.0, 0.0);
     Fsop->gemm(false, false, temp, Smhalf, 1.0, 0.0);

     // Obtain the orbitals
     Fsop->diagonalize(Cmop, e_orb, cutoff);
     CmoA->gemm(false, false, Smhalf, Cmop, 1.0, 0.0);
     if (reference_ == "UNRESTRICTED") CmoB->copy(CmoA);

     // memfree
     Fso.reset();
     Fsop.reset();
     Cmop.reset();
     temp.reset();
     Uso.reset();
     Smhalf.reset();
     e_orb.reset();
     DiagS.reset();

     // build mo coeff blocks
     mo_coeff_blocks();

}// end of gwh

//=======================================================
//          Canonic
//=======================================================
void DFOCC::canonic()
{
	SharedTensor2d UeigA = std::shared_ptr<Tensor2d>(new Tensor2d("UooA", nmo_, nmo_));
	SharedTensor1d eigA = std::shared_ptr<Tensor1d>(new Tensor1d("epsilon <A|A>", nmo_));

	// Diagonalize Fock
	FockA->diagonalize(UeigA, eigA, cutoff);

        // Build U
	UorbA->zero();
        UorbA->copy(UeigA);

        // Get new MOs
        SharedTensor2d Ca_new = std::shared_ptr<Tensor2d>(new Tensor2d("New alpha MO coefficients", nso_, nmo_));
	Ca_new->gemm(false, false, CmoA, UorbA, 1.0, 0.0);
	CmoA->copy(Ca_new);
	Ca_new.reset();

        // memfree
        UeigA.reset();
	eigA.reset();

//==========================================================================================
//========================= UHF REFERENCE ==================================================
//==========================================================================================
     if (reference_ == "UNRESTRICTED") {
       	SharedTensor2d UeigB = std::shared_ptr<Tensor2d>(new Tensor2d("UeigB", nmo_, nmo_));
	SharedTensor1d eigB = std::shared_ptr<Tensor1d>(new Tensor1d("epsilon <a|a>", nmo_));

	// Diagonalize Fock
	FockB->diagonalize(UeigB, eigB, cutoff);

        // Build U
	UorbB->zero();
        UorbB->copy(UeigB);

        // Get new MOs
        SharedTensor2d Cb_new = std::shared_ptr<Tensor2d>(new Tensor2d("New beta MO coefficients", nso_, nmo_));
	Cb_new->gemm(false, false, CmoB, UorbB, 1.0, 0.0);
	CmoB->copy(Cb_new);
	Cb_new.reset();

        // mem free
        UeigB.reset();
	eigB.reset();

     }// end uhf

     // build mo coeff blocks
     mo_coeff_blocks();

}// end of canonic

//=======================================================
//          Kappa
//=======================================================
void DFOCC::kappa_qchf()
{
//outfile->Printf("\n kappa_qchf is starting... \n");

    SharedTensor2d K, L;

if (reference_ == "RESTRICTED") {
    // Memalloc
    zvectorA = SharedTensor1d(new Tensor1d("Alpha Z-Vector", noccA * nvirA));
    zvec_newA = SharedTensor1d(new Tensor1d("Alpha New Z-Vector", noccA * nvirA));
    Minv_pcgA = SharedTensor1d(new Tensor1d("Alpha PCG M inverse", noccA * nvirA));
    sigma_pcgA = SharedTensor1d(new Tensor1d("Alpha PCG sigma", noccA * nvirA));
    r_pcgA = SharedTensor1d(new Tensor1d("Alpha PCG r", noccA * nvirA));
    r_pcg_newA = SharedTensor1d(new Tensor1d("Alpha PCG new r", noccA * nvirA));
    z_pcgA = SharedTensor1d(new Tensor1d("Alpha PCG z", noccA * nvirA));
    z_pcg_newA = SharedTensor1d(new Tensor1d("Alpha PCG new z", noccA * nvirA));
    p_pcgA = SharedTensor1d(new Tensor1d("Alpha PCG p", noccA * nvirA));
    p_pcg_newA = SharedTensor1d(new Tensor1d("Alpha PCG new p", noccA * nvirA));
    dr_pcgA = SharedTensor1d(new Tensor1d("Alpha PCG dr", noccA * nvirA));
    residualA = SharedTensor1d(new Tensor1d("Alpha Residual Vector", noccA * nvirA));

    // Build kappa0 and M
    for (int a = 0, ai = 0; a < nvirA; a++) {
         for (int i = 0; i < noccA; i++, ai++) {
              double value = FockA->get(a + noccA, a + noccA) - FockA->get(i,i);
              zvectorA->set(ai, -WorbA->get(a + noccA, i) / (2.0*value));
              Minv_pcgA->set(ai, 0.5/value);
         }
    }

    // Build S = A kappa_0
    sigma_rhf(sigma_pcgA, zvectorA);

    // Level Shift
    if (level_shift == "TRUE") sigma_pcgA->axpy(zvectorA, lshift_parameter);

    // Build r0
    for (int a = 0, ai = 0; a < nvirA; a++) {
         for (int i = 0; i < noccA; i++, ai++) {
              residualA->set(ai, -WorbA->get(a + noccA, i) - sigma_pcgA->get(ai));
         }
    }
    r_pcgA->copy(residualA);

    // Build z0
    z_pcgA->dirprd(Minv_pcgA, r_pcgA);

    // Build p0
    p_pcgA->copy(z_pcgA);

    // Call Orbital Response Solver
    orb_resp_pcg_rhf();

    // Memfree
    zvec_newA.reset();
    Minv_pcgA.reset();
    sigma_pcgA.reset();
    r_pcgA.reset();
    r_pcg_newA.reset();
    z_pcgA.reset();
    z_pcg_newA.reset();
    p_pcgA.reset();
    p_pcg_newA.reset();
    dr_pcgA.reset();
    residualA.reset();

    // Build kappa for VO block
    #pragma omp parallel for
    for (int x = 0; x < nidpA; x++) {
         int p = idprowA->get(x);
	 int q = idpcolA->get(x);
         if (p >= noccA && q < noccA) {
             int ai = vo_idxAA->get(p-noccA,q);
	     kappaA->set(x, zvectorA->get(ai));
         }
    }
    zvectorA.reset();

    // If LINEQ FAILED!
    if (pcg_conver == 0) {
        outfile->Printf("\tWarning!!! PCG did NOT converged in %2d iterations, switching to an approximately diagonal MO Hessian. \n", itr_pcg);

        // Build kappa again
        #pragma omp parallel for
        for (int x = 0; x < nidpA; x++) {
	    int p = idprowA->get(x);
	    int q = idpcolA->get(x);
            double value = FockA->get(p,p) - FockA->get(q,q);
	    kappaA->set(x, -wogA->get(x)/(2.0*value));
        }
    } // end if pcg_conver = 0

        // find biggest_kappa
	biggest_kappaA=0;
	for (int i=0; i<nidpA;i++) {
	    if (fabs(kappaA->get(i)) > biggest_kappaA) biggest_kappaA=fabs(kappaA->get(i));
	}

        // Scale
	if (biggest_kappaA > step_max) {
	    for (int i=0; i<nidpA;i++) kappaA->set(i, kappaA->get(i) *(step_max/biggest_kappaA));
	}

        // find biggest_kappa again
	if (biggest_kappaA > step_max)
	{
	  biggest_kappaA=0;
	  for (int i=0; i<nidpA;i++)
	  {
	      if (fabs(kappaA->get(i)) > biggest_kappaA)
	      {
		  biggest_kappaA = fabs(kappaA->get(i));
	      }
	  }
	}

        // norm
	rms_kappaA=0;
	rms_kappaA = kappaA->rms();

        // print
        if(print_ > 2) kappaA->print();

}// end if (reference_ == "RESTRICTED")

else if (reference_ == "UNRESTRICTED") {
    nidp_tot = nidpA + nidpB;

    // Memalloc
    zvector = SharedTensor1d(new Tensor1d("UHF Z-Vector", nidp_tot));
    zvec_new = SharedTensor1d(new Tensor1d("New UHF Z-Vector", nidp_tot));
    Minv_pcg = SharedTensor1d(new Tensor1d("PCG M inverse", nidp_tot));
    sigma_pcg = SharedTensor1d(new Tensor1d("PCG sigma", nidp_tot));
    r_pcg = SharedTensor1d(new Tensor1d("PCG r", nidp_tot));
    r_pcg_new = SharedTensor1d(new Tensor1d("PCG new r", nidp_tot));
    z_pcg = SharedTensor1d(new Tensor1d("PCG z", nidp_tot));
    z_pcg_new = SharedTensor1d(new Tensor1d("PCG new z", nidp_tot));
    p_pcg = SharedTensor1d(new Tensor1d("PCG p", nidp_tot));
    p_pcg_new = SharedTensor1d(new Tensor1d("PCG new p", nidp_tot));
    dr_pcg = SharedTensor1d(new Tensor1d("PCG dr", nidp_tot));
    residual = SharedTensor1d(new Tensor1d("Residual Vector", nidp_tot));
    zvectorA = SharedTensor1d(new Tensor1d("Alpha Z-Vector", noccA * nvirA));
    zvectorB = SharedTensor1d(new Tensor1d("Beta Z-Vector", noccB * nvirB));
    sigma_pcgA = SharedTensor1d(new Tensor1d("Alpha PCG sigma", noccA * nvirA));
    sigma_pcgB = SharedTensor1d(new Tensor1d("Beta PCG sigma", noccB * nvirB));
    p_pcgA = SharedTensor1d(new Tensor1d("Alpha PCG p", noccA * nvirA));
    p_pcgB = SharedTensor1d(new Tensor1d("Beta PCG p", noccB * nvirB));

    // Build kappa0 and M
    // alpha
    for (int a = 0, ai = 0; a < nvirA; a++) {
         for (int i = 0; i < noccA; i++, ai++) {
              double value = FockA->get(a + noccA, a + noccA) - FockA->get(i,i);
              zvectorA->set(ai, -WorbA->get(a + noccA, i) / (2.0*value));
              Minv_pcg->set(ai, 0.5/value);
         }
    }

    // beta
    for (int a = 0, ai = 0; a < nvirB; a++) {
         for (int i = 0; i < noccB; i++, ai++) {
              double value = FockB->get(a + noccB, a + noccB) - FockB->get(i,i);
              zvectorB->set(ai, -WorbB->get(a + noccB, i) / (2.0*value));
              Minv_pcg->set(ai + nidpA, 0.5/value);
         }
    }

    // Form initial zvector vector
    for (int ai = 0; ai < nidpA; ai++) zvector->set(ai, zvectorA->get(ai));
    for (int ai = 0; ai < nidpB; ai++) zvector->set(ai + nidpA, zvectorB->get(ai));

    // Build S = A kappa_0
    sigma_uhf(sigma_pcgA, sigma_pcgB, zvectorA, zvectorB);

    // Level Shift
    if (level_shift == "TRUE") {
        sigma_pcgA->axpy(zvectorA, lshift_parameter);
        sigma_pcgB->axpy(zvectorB, lshift_parameter);
    }

    // Form sigma vector
    for (int ai = 0; ai < nidpA; ai++) sigma_pcg->set(ai, sigma_pcgA->get(ai));
    for (int ai = 0; ai < nidpB; ai++) sigma_pcg->set(ai + nidpA, sigma_pcgB->get(ai));

    // Build r0
    // alpha
    for (int a = 0, ai = 0; a < nvirA; a++) {
         for (int i = 0; i < noccA; i++, ai++) {
              residual->set(ai, -WorbA->get(a + noccA, i));
         }
    }
    // beta
    for (int a = 0, ai = 0; a < nvirB; a++) {
         for (int i = 0; i < noccB; i++, ai++) {
              residual->set(ai + nidpA, -WorbB->get(a + noccB, i));
         }
    }
    residual->subtract(sigma_pcg);
    r_pcg->copy(residual);

    // Build z0
    z_pcg->dirprd(Minv_pcg, r_pcg);

    // Build p0
    p_pcg->copy(z_pcg);

    // Form initial pA and pB vectors
    for (int ai = 0; ai < nidpA; ai++) p_pcgA->set(ai, p_pcg->get(ai));
    for (int ai = 0; ai < nidpB; ai++) p_pcgB->set(ai, p_pcg->get(ai + nidpA));

    // Call Orbital Response Solver
    orb_resp_pcg_uhf();

    // Memfree alpha
    zvec_new.reset();
    Minv_pcg.reset();
    sigma_pcg.reset();
    sigma_pcgA.reset();
    sigma_pcgB.reset();
    r_pcg.reset();
    r_pcg_new.reset();
    z_pcg.reset();
    z_pcg_new.reset();
    p_pcg.reset();
    p_pcgA.reset();
    p_pcgB.reset();
    p_pcg_new.reset();
    dr_pcg.reset();
    residual.reset();

    // Build zvector for VO block
    for (int a = 0, ai = 0; a < nvirA; a++) {
         for (int i = 0; i < noccA; i++, ai++) {
	      zvectorA->set(ai, zvector->get(ai));
         }
    }

    // Build zvector for vo block
    for (int a = 0, ai = 0; a < nvirB; a++) {
         for (int i = 0; i < noccB; i++, ai++) {
	      zvectorB->set(ai, zvector->get(ai + nidpA));
         }
    }
    zvector.reset();

    // Build kappa for VO block
    // alpha
    #pragma omp parallel for
    for (int x = 0; x < nidpA; x++) {
         int p = idprowA->get(x);
	 int q = idpcolA->get(x);
         if (p >= noccA && q < noccA) {
             int ai = vo_idxAA->get(p-noccA,q);
	     kappaA->set(x, zvectorA->get(ai));
         }
    }
    zvectorA.reset();
    // beta
    #pragma omp parallel for
    for (int x = 0; x < nidpB; x++) {
         int p = idprowB->get(x);
	 int q = idpcolB->get(x);
         if (p >= noccB && q < noccB) {
             int ai = vo_idxBB->get(p-noccB,q);
	     kappaB->set(x, zvectorB->get(ai));
         }
    }
    zvectorB.reset();

    // If LINEQ FAILED!
    if (pcg_conver == 0) {
        outfile->Printf("\tWarning!!! PCG did NOT converged in %2d iterations, switching to an approximately diagonal MO Hessian. \n", itr_pcg);

	// alpha
        #pragma omp parallel for
        for(int x = 0; x < nidpA; x++) {
	    int p = idprowA->get(x);
	    int q = idpcolA->get(x);
            double value = FockA->get(p,p) - FockA->get(q,q);
	    kappaA->set(x, -wogA->get(x)/(2.0*value));
        }

	// beta
        #pragma omp parallel for
        for(int x = 0; x < nidpB; x++) {
	    int p = idprowB->get(x);
	    int q = idpcolB->get(x);
            double value = FockA->get(p,p) - FockA->get(q,q);
	    kappaB->set(x, -wogB->get(x)/(2.0*value));
        }
    } // end if pcg_conver = 0

        // find biggest_kappa
	biggest_kappaA=0;
	for (int i=0; i<nidpA;i++) {
	    if (fabs(kappaA->get(i)) > biggest_kappaA) biggest_kappaA=fabs(kappaA->get(i));
	}

	biggest_kappaB=0;
	for (int i=0; i<nidpB;i++){
	    if (fabs(kappaB->get(i)) > biggest_kappaB) biggest_kappaB=fabs(kappaB->get(i));
	}

        // Scale
	if (biggest_kappaA > step_max) {
	    for (int i=0; i<nidpA;i++) kappaA->set(i, kappaA->get(i) *(step_max/biggest_kappaA));
	}

	if (biggest_kappaB > step_max) {
	    for (int i=0; i<nidpB;i++) kappaB->set(i, kappaB->get(i) *(step_max/biggest_kappaB));
	}

        // find biggest_kappa again
	if (biggest_kappaA > step_max)
	{
	  biggest_kappaA=0;
	  for (int i=0; i<nidpA;i++)
	  {
	      if (fabs(kappaA->get(i)) > biggest_kappaA)
	      {
		  biggest_kappaA = fabs(kappaA->get(i));
	      }
	  }
	}

	if (biggest_kappaB > step_max)
	{
	  biggest_kappaB=0;
	  for (int i=0; i<nidpB;i++)
	  {
	      if (fabs(kappaB->get(i)) > biggest_kappaB)
	      {
		  biggest_kappaB=fabs(kappaB->get(i));
	      }
	  }
	}

        // norm
	rms_kappaA=0;
	rms_kappaB=0;
	rms_kappaA = kappaA->rms();
	rms_kappaB = kappaB->rms();

        // print
        if(print_ > 2){
          kappaA->print();
          kappaB->print();
        }
}// end if (reference_ == "UNRESTRICTED")
 //outfile->Printf("\n kappa_qchf done. \n");
}// end kappa_qchf

}} // End Namespaces
