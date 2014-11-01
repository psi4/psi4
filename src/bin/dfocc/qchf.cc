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

#include <libqt/qt.h>
#include "defines.h"
#include "dfocc.h"

using namespace psi;
using namespace std;


namespace psi{ namespace dfoccwave{
  
void DFOCC::qchf()
{
   
fprintf(outfile,"\n");      
fprintf(outfile," ============================================================================== \n");    
fprintf(outfile," ================ Performing QCHF iterations... =============================== \n");  
fprintf(outfile," ============================================================================== \n");
fprintf(outfile, "\t            QCHF \n");
fprintf(outfile, "\t   ------------------------------ \n");
//fprintf(outfile, " Iter       E_total           DE           RMS MO Grad \n");
//fprintf(outfile, " ----    ---------------    ----------     ----------- \n");
fprintf(outfile, " Iter       E_total           DE           RMS MO Grad      MAX MO Grad  \n");
fprintf(outfile, " ----    ---------------    ----------     -----------      -----------  \n");
fflush(outfile);

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
        kappa_orb_resp_pcg();
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
	
fprintf(outfile," %3d     %12.10f  %12.2e   %12.2e     %12.2e \n",itr_occ,Eref,DE,rms_wog,biggest_mograd);
fflush(outfile);

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
fprintf(outfile,"\n");
fprintf(outfile," ============================================================================== \n");
fprintf(outfile," ======================== QCHF ITERATIONS ARE CONVERGED ======================= \n");
fprintf(outfile," ============================================================================== \n");
fflush(outfile);

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
  fprintf(outfile,"\n\tQCHF IS NOT CONVERGED IN %2d ITERATIONS ========== \n", mo_maxiter);
  fflush(outfile);
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
            fprintf(outfile,"\tThere is not any non-redundant orbital rotation pair! \n");
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
        fprintf(outfile,"\tThere is not any non-redundant orbital rotation pair! \n");
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
     SharedTensor1d e_orb = boost::shared_ptr<Tensor1d>(new Tensor1d("epsilon <n|n>", nso_));
     SharedTensor1d DiagS = boost::shared_ptr<Tensor1d>(new Tensor1d("Diag S", nso_));

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

if (reference_ == "RESTRICTED") { 
        // Build Cocc
        for (int mu = 0; mu < nso_; mu++) {
             for (int i = 0; i < noccA; i++) {
                 CoccA->set(mu, i, CmoA->get(mu, i));
             }
        }

        // Build Cvir
        for (int mu = 0; mu < nso_; mu++) {
             for (int a = 0; a < nvirA; a++) {
                 CvirA->set(mu, a, CmoA->get(mu, a + noccA));
             }
        }
 
        // Build active Caocc
        for (int mu = 0; mu < nso_; mu++) {
             for (int i = 0; i < naoccA; i++) {
                 CaoccA->set(mu, i, CmoA->get(mu, i + nfrzc));
             }
        }

        // Build active Cvir
        for (int mu = 0; mu < nso_; mu++) {
             for (int a = 0; a < navirA; a++) {
                 CavirA->set(mu, a, CmoA->get(mu, a + noccA));
             }
        }
}// end if (reference_ == "RESTRICTED") 

else if (reference_ == "UNRESTRICTED") {
        // Build Cocc
        // alpha
        for (int mu = 0; mu < nso_; mu++) {
             for (int i = 0; i < noccA; i++) {
                 CoccA->set(mu, i, CmoA->get(mu, i));
             }
        }

        // beta
        for (int mu = 0; mu < nso_; mu++) {
             for (int i = 0; i < noccB; i++) {
                 CoccB->set(mu, i, CmoB->get(mu, i));
             }
        }

        // Build Cvir
        // alpha
        for (int mu = 0; mu < nso_; mu++) {
             for (int a = 0; a < nvirA; a++) {
                 CvirA->set(mu, a, CmoA->get(mu, a + noccA));
             }
        }
 
        // beta
        for (int mu = 0; mu < nso_; mu++) {
             for (int a = 0; a < nvirB; a++) {
                 CvirB->set(mu, a, CmoB->get(mu, a + noccB));
             }
        }

        // Build active Caocc
        // alpha
        for (int mu = 0; mu < nso_; mu++) {
             for (int i = 0; i < naoccA; i++) {
                 CaoccA->set(mu, i, CmoA->get(mu, i + nfrzc));
             }
        }

        // beta
        for (int mu = 0; mu < nso_; mu++) {
             for (int i = 0; i < naoccB; i++) {
                 CaoccB->set(mu, i, CmoB->get(mu, i + nfrzc));
             }
        }

        // Build active Cvir
        // alpha
        for (int mu = 0; mu < nso_; mu++) {
             for (int a = 0; a < navirA; a++) {
                 CavirA->set(mu, a, CmoA->get(mu, a + noccA));
             }
        }
 
        // beta
        for (int mu = 0; mu < nso_; mu++) {
             for (int a = 0; a < navirB; a++) {
                 CavirB->set(mu, a, CmoB->get(mu, a + noccB));
             }
        }

}// end if (reference_ == "UNRESTRICTED") 

}// end of gwh

//=======================================================
//          Canonic
//=======================================================          
void DFOCC::canonic()
{
	SharedTensor2d UeigA = boost::shared_ptr<Tensor2d>(new Tensor2d("UooA", nmo_, nmo_));
	SharedTensor1d eigA = boost::shared_ptr<Tensor1d>(new Tensor1d("epsilon <A|A>", nmo_));

	// Diagonalize Fock  
	FockA->diagonalize(UeigA, eigA, cutoff);

        // Build U	
	UorbA->zero();
        UorbA->copy(UeigA);

        // Get new MOs
        SharedTensor2d Ca_new = boost::shared_ptr<Tensor2d>(new Tensor2d("New alpha MO coefficients", nso_, nmo_));
	Ca_new->gemm(false, false, CmoA, UorbA, 1.0, 0.0); 
	CmoA->copy(Ca_new);
	Ca_new.reset();

        // memfree
        UeigA.reset();
	eigA.reset();

        // Build Cocc
        for (int mu = 0; mu < nso_; mu++) {
             for (int i = 0; i < noccA; i++) {
                 CoccA->set(mu, i, CmoA->get(mu, i));
             }
        }

        // Build Cvir
        for (int mu = 0; mu < nso_; mu++) {
             for (int a = 0; a < nvirA; a++) {
                 CvirA->set(mu, a, CmoA->get(mu, a + noccA));
             }
        }
 
        // Build active Caocc
        for (int mu = 0; mu < nso_; mu++) {
             for (int i = 0; i < naoccA; i++) {
                 CaoccA->set(mu, i, CmoA->get(mu, i + nfrzc));
             }
        }

        // Build active Cvir
        for (int mu = 0; mu < nso_; mu++) {
             for (int a = 0; a < navirA; a++) {
                 CavirA->set(mu, a, CmoA->get(mu, a + noccA));
             }
        }
//==========================================================================================
//========================= UHF REFERENCE ==================================================
//==========================================================================================
     if (reference_ == "UNRESTRICTED") {
       	SharedTensor2d UeigB = boost::shared_ptr<Tensor2d>(new Tensor2d("UeigB", nmo_, nmo_));
	SharedTensor1d eigB = boost::shared_ptr<Tensor1d>(new Tensor1d("epsilon <a|a>", nmo_));

	// Diagonalize Fock  
	FockB->diagonalize(UeigB, eigB, cutoff);

        // Build U	
	UorbB->zero();
        UorbB->copy(UeigB);

        // Get new MOs
        SharedTensor2d Cb_new = boost::shared_ptr<Tensor2d>(new Tensor2d("New beta MO coefficients", nso_, nmo_));
	Cb_new->gemm(false, false, CmoB, UorbB, 1.0, 0.0); 
	CmoB->copy(Cb_new);
	Cb_new.reset();

        // mem free
        UeigB.reset();
	eigB.reset();

        // Build Cocc
        for (int mu = 0; mu < nso_; mu++) {
             for (int i = 0; i < noccB; i++) {
                 CoccB->set(mu, i, CmoB->get(mu, i));
             }
        }

        // Build Cvir
        for (int mu = 0; mu < nso_; mu++) {
             for (int a = 0; a < nvirB; a++) {
                 CvirB->set(mu, a, CmoB->get(mu, a + noccB));
             }
        }
 
        // Build active Caocc
        for (int mu = 0; mu < nso_; mu++) {
             for (int i = 0; i < naoccB; i++) {
                 CaoccB->set(mu, i, CmoB->get(mu, i + nfrzc));
             }
        }

        // Build active Cvir
        for (int mu = 0; mu < nso_; mu++) {
             for (int a = 0; a < navirB; a++) {
                 CavirB->set(mu, a, CmoB->get(mu, a + noccB));
             }
        }
     }// end uhf	

}// end of canonic

}} // End Namespaces


