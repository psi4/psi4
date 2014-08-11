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

/** Standard library includes */
#include <fstream>
#include <psifiles.h>
#include <libiwl/iwl.hpp>
#include <libqt/qt.h>
#include "dfocc.h"

using namespace boost;
using namespace psi;
using namespace std;

namespace psi{ namespace dfoccwave{

void DFOCC::get_moinfo()
{      
  //outfile->Printf("\n get_moinfo is starting... \n"); 
//===========================================================================================
//========================= RHF =============================================================
//===========================================================================================
if (reference_ == "RESTRICTED") {
	// Read in mo info
        nso_     = reference_wavefunction_->nso();
        nirrep_  = reference_wavefunction_->nirrep();
        nmo_     = reference_wavefunction_->nmo();
        nmopi_   = reference_wavefunction_->nmopi();
        nsopi_   = reference_wavefunction_->nsopi();
        doccpi_  = reference_wavefunction_->doccpi();
        frzcpi_  = reference_wavefunction_->frzcpi();
        frzvpi_  = reference_wavefunction_->frzvpi();
        natom   = molecule_->natom();

	// Read in nuclear repulsion energy
	Enuc = Process::environment.molecule()->nuclear_repulsion_energy();
	
	// Read SCF energy
        Escf=reference_wavefunction_->reference_energy();
	Eref=Escf;
	Eelec=Escf-Enuc;
	
        // figure out number of MO spaces
        nfrzc = frzcpi_[0];
        nfrzv = frzvpi_[0];
	noccA = doccpi_[0];
	nvirA = nmo_ - noccA;   	      // Number of virtual orbitals
	naoccA = noccA - nfrzc; 	      // Number of active occupied orbitals
	navirA = nvirA - nfrzv;               // Number of active virtual orbitals
        nocc2AA = noccA * noccA;              // Number of all OCC-OCC pairs 
        nvir2AA = nvirA * nvirA;              // Number of all VIR-VIR pairs 
        naocc2AA = naoccA * naoccA;           // Number of active OCC-OCC pairs 
        navir2AA = navirA * navirA;           // Number of active VIR-VIR pairs 
        navoAA = naoccA * navirA;             // Number of active OCC-VIR pairs
        nvoAA = noccA * nvirA;                // Number of all OCC-VIR pairs
	namo = nmo_- nfrzc - nfrzv; 	      // Number of active  orbitals
	npop = nmo_ - nfrzv;         	      // Number of populated orbitals
	ntri_so = 0.5*nso_*(nso_+1);
        ntri = 0.5*nmo_*(nmo_+1);
	dimtei = 0.5*ntri*(ntri+1);
        nso2_ = nso_ * nso_;
    
/********************************************************************************************/
/************************** Read orbital coefficients ***************************************/
/********************************************************************************************/
        // Read orbital energies
        epsilon_a_ = reference_wavefunction_->epsilon_a();

        // Build Initial fock matrix
	FockA = SharedTensor2d(new Tensor2d("MO-basis alpha Fock matrix", nmo_, nmo_));
        for(int i = 0; i < noccA; ++i) FockA->set(i, i, epsilon_a_->get(0, i));
        for(int a = 0; a < nvirA; ++a) FockA->set(a + noccA, a + noccA, epsilon_a_->get(0, a + noccA));
        if (print_ > 2) FockA->print();

        // Fock Blocks
        FooA = SharedTensor2d(new Tensor2d("Fock <O|O>", noccA, noccA));
        FovA = SharedTensor2d(new Tensor2d("Fock <O|V>", noccA, nvirA));
        FvoA = SharedTensor2d(new Tensor2d("Fock <V|O>", nvirA, noccA));
        FvvA = SharedTensor2d(new Tensor2d("Fock <V|V>", nvirA, nvirA));
        FooA->form_oo(FockA);
        FvoA->form_vo(FockA);
        FovA = FvoA->transpose();
        FvvA->form_vv(noccA, FockA);

        // Figure out OO Scaling factor
        if (options_["OO_SCALE"].has_changed()) {
            msd_oo_scale=options_.get_double("OO_SCALE");
        }
        else {
            //msd_oo_scale = ( FockA->get(noccA,noccA) - FockA->get(noccA-1,noccA-1) ) / ( FockA->get(1,1) - FockA->get(0,0) ); 
            msd_oo_scale = ( FockA->get(noccA,noccA) - FockA->get(noccA-1,noccA-1) ) / ( FockA->get(nfrzc+1,nfrzc+1) - FockA->get(nfrzc,nfrzc) ); 
            msd_oo_scale /= 3.0;
            if (hess_type == "APPROX_DIAG_HF" || hess_type == "APPROX_DIAG_EKT") {
                outfile->Printf("\tOO Scale is changed to: %12.10f\n", msd_oo_scale);
                
            }
        }

        // Read orbital coefficients from chkpt
	Ca_ = SharedMatrix(reference_wavefunction_->Ca());
        CmoA = SharedTensor2d(new Tensor2d("Alpha MO Coefficients", nso_, nmo_));
        CmoA->set(Ca_);
        if (orb_opt_ == "TRUE") {
            Cmo_refA = SharedTensor2d(new Tensor2d("Alpha Reference MO Coefficients", nso_, nmo_));
            Cmo_refA->copy(CmoA);
        }
        if (print_ > 2) CmoA->print();

        // Build Cocc
        CoccA = SharedTensor2d(new Tensor2d("Alpha C(mu,i)", nso_, noccA));
        for (int mu = 0; mu < nso_; mu++) {
             for (int i = 0; i < noccA; i++) {
                 CoccA->set(mu, i, CmoA->get(mu, i));
             }
        }
        if (print_ > 2) CoccA->print();

        // Build Cvir
        CvirA = SharedTensor2d(new Tensor2d("Alpha C(mu,a)", nso_, nvirA));
        for (int mu = 0; mu < nso_; mu++) {
             for (int a = 0; a < nvirA; a++) {
                 CvirA->set(mu, a, CmoA->get(mu, a + noccA));
             }
        }
        if (print_ > 2) CvirA->print();
 

        // Build active Caocc
        CaoccA = SharedTensor2d(new Tensor2d("Alpha Active C(mu,i)", nso_, naoccA));
        for (int mu = 0; mu < nso_; mu++) {
             for (int i = 0; i < naoccA; i++) {
                 CaoccA->set(mu, i, CmoA->get(mu, i + nfrzc));
             }
        }
        if (print_ > 2) CaoccA->print();

        // Build active Cvir
        CavirA = SharedTensor2d(new Tensor2d("Alpha Active C(mu,a)", nso_, navirA));
        for (int mu = 0; mu < nso_; mu++) {
             for (int a = 0; a < navirA; a++) {
                 CavirA->set(mu, a, CmoA->get(mu, a + noccA));
             }
        }
        if (print_ > 2) CavirA->print();
}// if (reference_ == "RESTRICTED") 

//===========================================================================================
//========================= UHF =============================================================
//===========================================================================================
else if (reference_ == "UNRESTRICTED") {
	// Read in mo info
        nso_     = reference_wavefunction_->nso();
        nirrep_  = reference_wavefunction_->nirrep();
        nmo_     = reference_wavefunction_->nmo();
        nmopi_   = reference_wavefunction_->nmopi();
        nsopi_   = reference_wavefunction_->nsopi();
        doccpi_  = reference_wavefunction_->doccpi();
        soccpi_  = reference_wavefunction_->soccpi();
        frzcpi_  = reference_wavefunction_->frzcpi();
        frzvpi_  = reference_wavefunction_->frzvpi();
        natom   = molecule_->natom();

	// Read in nuclear repulsion energy
	Enuc = Process::environment.molecule()->nuclear_repulsion_energy();
	
	// Read SCF energy
        Escf=reference_wavefunction_->reference_energy();
	Eref=Escf;
	Eelec=Escf-Enuc;

        // figure out number of MO spaces
        nfrzc = frzcpi_[0];
        nfrzv = frzvpi_[0];
	noccA = doccpi_[0] + soccpi_[0];
	noccB = doccpi_[0];
	nvirA = nmo_ - noccA;   		// Number of virtual orbitals
	nvirB = nmo_ - noccB;   	        // Number of virtual orbitals
	naoccA = noccA - nfrzc; 		// Number of active occupied orbitals
	naoccB = noccB - nfrzc; 	        // Number of active occupied orbitals
	navirA = nvirA - nfrzv; 		// Number of active virtual orbitals
	navirB = nvirB - nfrzv; 	        // Number of active virtual orbitals
        nocc2AA = noccA * noccA;                // Number of all OCC-OCC pairs 
        nocc2AB = noccA * noccB;                // Number of all OCC-OCC pairs 
        nocc2BB = noccB * noccB;                // Number of all OCC-OCC pairs 
        nvir2AA = nvirA * nvirA;                // Number of all VIR-VIR pairs 
        nvir2AB = nvirA * nvirB;                // Number of all VIR-VIR pairs 
        nvir2BB = nvirB * nvirB;                // Number of all VIR-VIR pairs 
        naocc2AA = naoccA * naoccA;              // Number of active OCC-OCC pairs 
        naocc2AB = naoccA * naoccB;              // Number of active OCC-OCC pairs 
        naocc2BB = naoccB * naoccB;              // Number of active OCC-OCC pairs 
        navir2AA = navirA * navirA;              // Number of active VIR-VIR pairs 
        navir2AB = navirA * navirB;              // Number of active VIR-VIR pairs 
        navir2BB = navirB * navirB;              // Number of active VIR-VIR pairs 
        navoAA = naoccA * navirA;                // Number of active OCC-VIR pairs
        navoBA = naoccA * navirB;                // Number of active OCC-VIR pairs
        navoAB = naoccB * navirA;                // Number of active OCC-VIR pairs
        navoBB = naoccB * navirB;                // Number of active OCC-VIR pairs
        nvoAA = noccA * nvirA;                   // Number of all OCC-VIR pairs
        nvoBA = noccA * nvirB;                   // Number of all OCC-VIR pairs
        nvoAB = noccB * nvirA;                   // Number of all OCC-VIR pairs
        nvoBB = noccB * nvirB;                   // Number of all OCC-VIR pairs
	namo = nmo_- nfrzc - nfrzv; 	         // Number of active  orbitals
	npop = nmo_ - nfrzv;         	         // Number of populated orbitals
	ntri_so = 0.5*nso_*(nso_+1);
        ntri = 0.5*nmo_*(nmo_+1);
	dimtei = 0.5*ntri*(ntri+1);
        nso2_ = nso_ * nso_;

/********************************************************************************************/
/************************** Read orbital coefficients ***************************************/
/********************************************************************************************/
        // Read orbital energies
        epsilon_a_ = reference_wavefunction_->epsilon_a();
        epsilon_b_ = reference_wavefunction_->epsilon_b();

        // Build Initial fock matrix
	FockA = SharedTensor2d(new Tensor2d("MO-basis alpha Fock matrix", nmo_, nmo_));
        for(int i = 0; i < noccA; ++i) FockA->set(i, i, epsilon_a_->get(0, i));
        for(int a = 0; a < nvirA; ++a) FockA->set(a + noccA, a + noccA, epsilon_a_->get(0, a + noccA));
        if (print_ > 2) FockA->print();

	FockB = SharedTensor2d(new Tensor2d("MO-basis beta Fock matrix", nmo_, nmo_));
        for(int i = 0; i < noccB; ++i) FockB->set(i, i, epsilon_b_->get(0, i));
        for(int a = 0; a < nvirB; ++a) FockB->set(a + noccB, a + noccB, epsilon_b_->get(0, a + noccB));
        if (print_ > 2) FockB->print();

        // Fock Blocks
        FooA = SharedTensor2d(new Tensor2d("Fock <O|O>", noccA, noccA));
        FooB = SharedTensor2d(new Tensor2d("Fock <o|o>", noccB, noccB));
        FovA = SharedTensor2d(new Tensor2d("Fock <O|V>", noccA, nvirA));
        FovB = SharedTensor2d(new Tensor2d("Fock <o|v>", noccB, nvirB));
        FvoA = SharedTensor2d(new Tensor2d("Fock <V|O>", nvirA, noccA));
        FvoB = SharedTensor2d(new Tensor2d("Fock <v|o>", nvirB, noccB));
        FvvA = SharedTensor2d(new Tensor2d("Fock <V|V>", nvirA, nvirA));
        FvvB = SharedTensor2d(new Tensor2d("Fock <v|v>", nvirB, nvirB));
        FooA->form_oo(FockA);
        FooB->form_oo(FockB);
        FvoA->form_vo(FockA);
        FvoB->form_vo(FockB);
        FovA = FvoA->transpose();
        FovB = FvoB->transpose();
        FvvA->form_vv(noccA, FockA);
        FvvB->form_vv(noccB, FockB);

        // Figure out OO Scaling factor
        if (options_["OO_SCALE"].has_changed()) {
            msd_oo_scale=options_.get_double("OO_SCALE");
        }
        else {
            //double scaleA = ( FockA->get(noccA,noccA) - FockA->get(noccA-1,noccA-1) ) / ( FockA->get(1,1) - FockA->get(0,0) ); 
            //double scaleB = ( FockB->get(noccB,noccB) - FockB->get(noccB-1,noccB-1) ) / ( FockB->get(1,1) - FockB->get(0,0) ); 
            double scaleA = ( FockA->get(noccA,noccA) - FockA->get(noccA-1,noccA-1) ) / ( FockA->get(nfrzc+1,nfrzc+1) - FockA->get(nfrzc,nfrzc) ); 
            double scaleB = ( FockB->get(noccB,noccB) - FockB->get(noccB-1,noccB-1) ) / ( FockB->get(nfrzc+1,nfrzc+1) - FockB->get(nfrzc,nfrzc) ); 
            scaleA /= 3.0;
            scaleB /= 3.0;
            msd_oo_scale = 0.5 * (scaleA + scaleB);
            if (hess_type == "APPROX_DIAG_HF" || hess_type == "APPROX_DIAG_EKT") {
                outfile->Printf("\tOO Scale is changed to: %12.10f\n", msd_oo_scale);
                
            }
        }

        // Read orbital coefficients from chkpt
	Ca_ = SharedMatrix(reference_wavefunction_->Ca());
        CmoA = SharedTensor2d(new Tensor2d("Alpha MO Coefficients", nso_, nmo_));
        CmoA->set(Ca_);
        if (print_ > 2) CmoA->print();

        Cb_ = SharedMatrix(reference_wavefunction_->Cb());
        CmoB = SharedTensor2d(new Tensor2d("Beta MO Coefficients", nso_, nmo_));
        CmoB->set(Cb_);
        if (orb_opt_ == "TRUE") {
            Cmo_refA = SharedTensor2d(new Tensor2d("Alpha Reference MO Coefficients", nso_, nmo_));
            Cmo_refB = SharedTensor2d(new Tensor2d("Beta Reference MO Coefficients", nso_, nmo_));
            Cmo_refA->copy(CmoA);
            Cmo_refB->copy(CmoB);
        }
        if (print_ > 2) CmoB->print();

        // Build Cocc
        CoccA = SharedTensor2d(new Tensor2d("Alpha C(mu,i)", nso_, noccA));
        for (int mu = 0; mu < nso_; mu++) {
             for (int i = 0; i < noccA; i++) {
                 CoccA->set(mu, i, CmoA->get(mu, i));
             }
        }
        if (print_ > 2) CoccA->print();

        CoccB = SharedTensor2d(new Tensor2d("Beta C(mu,i)", nso_, noccB));
        for (int mu = 0; mu < nso_; mu++) {
             for (int i = 0; i < noccB; i++) {
                 CoccB->set(mu, i, CmoB->get(mu, i));
             }
        }
        if (print_ > 2) CoccB->print();

        // Build Cvir
        CvirA = SharedTensor2d(new Tensor2d("Alpha C(mu,a)", nso_, nvirA));
        for (int mu = 0; mu < nso_; mu++) {
             for (int a = 0; a < nvirA; a++) {
                 CvirA->set(mu, a, CmoA->get(mu, a + noccA));
             }
        }
        if (print_ > 2) CvirA->print();
 
        CvirB = SharedTensor2d(new Tensor2d("Beta C(mu,a)", nso_, nvirB));
        for (int mu = 0; mu < nso_; mu++) {
             for (int a = 0; a < nvirB; a++) {
                 CvirB->set(mu, a, CmoB->get(mu, a + noccB));
             }
        }
        if (print_ > 2) CvirB->print();

        // Build active Caocc
        CaoccA = SharedTensor2d(new Tensor2d("Alpha Active C(mu,i)", nso_, naoccA));
        for (int mu = 0; mu < nso_; mu++) {
             for (int i = 0; i < naoccA; i++) {
                 CaoccA->set(mu, i, CmoA->get(mu, i + nfrzc));
             }
        }
        if (print_ > 2) CaoccA->print();

        CaoccB = SharedTensor2d(new Tensor2d("Beta Active C(mu,i)", nso_, naoccB));
        for (int mu = 0; mu < nso_; mu++) {
             for (int i = 0; i < naoccB; i++) {
                 CaoccB->set(mu, i, CmoB->get(mu, i + nfrzc));
             }
        }
        if (print_ > 2) CaoccB->print();

        // Build active Cvir
        CavirA = SharedTensor2d(new Tensor2d("Alpha Active C(mu,a)", nso_, navirA));
        for (int mu = 0; mu < nso_; mu++) {
             for (int a = 0; a < navirA; a++) {
                 CavirA->set(mu, a, CmoA->get(mu, a + noccA));
             }
        }
        if (print_ > 2) CavirA->print();
 
        CavirB = SharedTensor2d(new Tensor2d("Beta Active C(mu,a)", nso_, navirB));
        for (int mu = 0; mu < nso_; mu++) {
             for (int a = 0; a < navirB; a++) {
                 CavirB->set(mu, a, CmoB->get(mu, a + noccB));
             }
        }
        if (print_ > 2) CavirB->print();

    if (reference == "ROHF") {
        Fa_ = SharedMatrix(reference_wavefunction_->Fa());
        Fb_ = SharedMatrix(reference_wavefunction_->Fb());
        FsoA = SharedTensor2d(new Tensor2d("SO-basis Alpha Fock Matrix", nso_, nso_));
        FsoB = SharedTensor2d(new Tensor2d("SO-basis Beta Fock Matrix", nso_, nso_));
        FsoA->set(Fa_);
        FsoB->set(Fb_);
        FockA->transform(FsoA, CmoA);
        FockB->transform(FsoB, CmoB);
    }
    
}// end else if (reference_ == "UNRESTRICTED") 

/********************************************************************************************/
/************************** Create all required matrice *************************************/
/********************************************************************************************/
        // Build Hso
	Hso_ = boost::shared_ptr<Matrix>(new Matrix("SO-basis One-electron Ints", nso_, nso_));
	Tso_ = boost::shared_ptr<Matrix>(new Matrix("SO-basis Kinetic Energy Ints", nso_, nso_));
	Vso_ = boost::shared_ptr<Matrix>(new Matrix("SO-basis Potential Energy Ints", nso_, nso_));
	Sso_ = boost::shared_ptr<Matrix>(new Matrix("SO-basis Overlap Ints", nso_, nso_));
	Hso_->zero();
	Tso_->zero();
	Vso_->zero();
	Sso_->zero();
	
	// Read SO-basis one-electron integrals
	double *so_ints = init_array(ntri_so);
        IWL::read_one(psio_.get(), PSIF_OEI, PSIF_SO_T, so_ints, ntri_so, 0, 0, "outfile");
        Tso_->set(so_ints);
        IWL::read_one(psio_.get(), PSIF_OEI, PSIF_SO_V, so_ints, ntri_so, 0, 0, "outfile");
        Vso_->set(so_ints);
        IWL::read_one(psio_.get(), PSIF_OEI, PSIF_SO_S, so_ints, ntri_so, 0, 0, "outfile");
        Sso_->set(so_ints);
        free(so_ints);
	Hso_->copy(Tso_); 
	Hso_->add(Vso_);
        Tso_.reset();
        Vso_.reset();
	Hso = SharedTensor2d(new Tensor2d("SO-basis One-electron Ints", nso_, nso_));
        Hso->set(Hso_);
        Hso_.reset();
	Sso = SharedTensor2d(new Tensor2d("SO-basis Overlap Ints", nso_, nso_));
        Sso->set(Sso_);
        Sso_.reset();

//outfile->Printf("\n get_moinfo is done. \n"); 
}// end get_moinfo
}} // End Namespaces




