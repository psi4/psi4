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

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <algorithm>
#include <vector>
#include <utility>

#include "psi4/libfunctional/superfunctional.h"
#include "psi4/libciomr/libciomr.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/libiwl/iwl.hpp"
#include "psi4/libqt/qt.h"

#include "psi4/libfock/jk.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/factory.h"
#include "cuhf.h"

using namespace std;
using namespace psi;


namespace psi { namespace scf {

CUHF::CUHF(SharedWavefunction ref_wfn, std::shared_ptr<SuperFunctional> func)
    : HF(ref_wfn, func, Process::environment.options, PSIO::shared_object())
{
    common_init();
}

CUHF::CUHF(SharedWavefunction ref_wfn, std::shared_ptr<SuperFunctional> func,
           Options& options, std::shared_ptr<PSIO> psio)
    : HF(ref_wfn, func, options, psio)
{
    common_init();
}

CUHF::~CUHF()
{
}

void CUHF::common_init()
{

    Drms_ = 0.0;

    Fa_         = SharedMatrix(factory_->create_matrix("F alpha"));
    Fb_         = SharedMatrix(factory_->create_matrix("F beta"));
    Fp_         = SharedMatrix(factory_->create_matrix("F charge"));
    Fm_         = SharedMatrix(factory_->create_matrix("F spin"));
    Da_         = SharedMatrix(factory_->create_matrix("SCF alpha density"));
    Db_         = SharedMatrix(factory_->create_matrix("SCF beta density"));
    Dp_         = SharedMatrix(factory_->create_matrix("D charge"));
    Dt_         = SharedMatrix(factory_->create_matrix("D total"));
    Da_old_     = SharedMatrix(factory_->create_matrix("D alpha old"));
    Db_old_     = SharedMatrix(factory_->create_matrix("D beta  old"));
    Dt_old_     = SharedMatrix(factory_->create_matrix("D total old"));
    Lagrangian_ = SharedMatrix(factory_->create_matrix("Lagrangian"));
    Ca_         = SharedMatrix(factory_->create_matrix("C alpha"));
    Cb_         = SharedMatrix(factory_->create_matrix("C beta"));
    Cno_        = SharedMatrix(factory_->create_matrix("C NOs"));
    Cno_temp_   = SharedMatrix(factory_->create_matrix("C NO temp"));
    J_          = SharedMatrix(factory_->create_matrix("J total"));
    Ka_         = SharedMatrix(factory_->create_matrix("K alpha"));
    Kb_         = SharedMatrix(factory_->create_matrix("K beta"));

    epsilon_a_ = SharedVector(factory_->create_vector());
    epsilon_b_ = SharedVector(factory_->create_vector());
    No_ = SharedVector(factory_->create_vector());
    same_a_b_dens_ = false;
    same_a_b_orbs_ = false;

    if (functional_->needs_xc()){
        throw PSIEXCEPTION("CUHF: Cannot compute XC components!");
    }
}

void CUHF::damp_update()
{
  Da_->scale(1.0 - damping_percentage_);
  Da_->axpy(damping_percentage_, Da_old_);
  Db_->scale(1.0 - damping_percentage_);
  Db_->axpy(damping_percentage_, Db_old_);
  Dt_->copy(Da_);
  Dt_->add(Db_);
}

void CUHF::finalize()
{
    // Form lagrangian
    for (int h=0; h<nirrep_; ++h) {
        for (int m=0; m<Lagrangian_->rowdim(h); ++m) {
            for (int n=0; n<Lagrangian_->coldim(h); ++n) {
                double sum = 0.0;
                for (int i=0; i<doccpi_[h]; ++i) {
                    sum += epsilon_a_->get(h, i) * Ca_->get(h, m, i) * Ca_->get(h, n, i)
                        +  epsilon_b_->get(h, i) * Cb_->get(h, m, i) * Cb_->get(h, n, i);
                }
                for (int i=doccpi_[h]; i<doccpi_[h]+soccpi_[h]; ++i)
                    sum += epsilon_a_->get(h, i) * Ca_->get(h, m, i) * Ca_->get(h, n, i);

                Lagrangian_->set(h, m, n, sum);
            }
        }
    }

    Dt_.reset();
    Da_old_.reset();
    Db_old_.reset();
    Dt_old_.reset();
    Dp_.reset();
    Fp_.reset();
    Fm_.reset();
    Cno_.reset();
    Cno_temp_.reset();
    No_.reset();

    HF::finalize();
}

void CUHF::save_density_and_energy()
{
    Da_old_->copy(Dt_);
    Db_old_->copy(Dt_);
    Dt_old_->copy(Dt_);
    Eold_ = E_;
}

void CUHF::form_G()
{
    // Push the C matrix on
    std::vector<SharedMatrix> & C = jk_->C_left();
    C.clear();
    C.push_back(Ca_subset("SO", "OCC"));
    C.push_back(Cb_subset("SO", "OCC"));

    // Run the JK object
    jk_->compute();

    // Pull the J and K matrices off
    const std::vector<SharedMatrix> & J = jk_->J();
    const std::vector<SharedMatrix> & K = jk_->K();
    J_->copy(J[0]);
    J_->add(J[1]);
    Ka_ = K[0];
    Kb_ = K[1];
}

void CUHF::save_information()
{
}

void CUHF::compute_spin_contamination()
{
    double dN = 0.0;

    for (int h =0; h < S_->nirrep(); h++) {
        int nbf = S_->colspi()[h];
        int nmo = Ca_->colspi()[h];
        int na = nalphapi_[h];
        int nb = nbetapi_[h];
        if (na == 0 || nb == 0 || nbf == 0 || nmo == 0)
            continue;

        SharedMatrix Ht (new Matrix("H Temp", nbf, nb));
        SharedMatrix Ft (new Matrix("F Temp", na, nb));

        double** Sp = S_->pointer(h);
        double** Cap = Ca_->pointer(h);
        double** Cbp = Cb_->pointer(h);
        double** Htp = Ht->pointer(0);
        double** Ftp = Ft->pointer(0);

        C_DGEMM('N','N',nbf,nb,nbf,1.0,Sp[0],nbf,Cbp[0],nmo,0.0,Htp[0],nb);
        C_DGEMM('T','N',na,nb,nbf,1.0,Cap[0],nmo,Htp[0],nb,0.0,Ftp[0],nb);

        for (long int ab = 0; ab < (long int)na*nb; ab++)
            dN += Ftp[0][ab]*Ftp[0][ab];
    }

    double dS = (double)nbeta_ - (double)dN;

    double nm = (nalpha_ - nbeta_) / 2.0;
    double S2 = nm * (nm + 1.0);

    outfile->Printf( "\n  @Spin Contamination Metric: %8.5F\n", dS);
      outfile->Printf( "  @S^2 Expected:              %8.5F\n", S2);
      outfile->Printf( "  @S^2 Observed:              %8.5F\n", S2 + dS);

}

bool CUHF::test_convergency()
{
    double ediff = E_ - Eold_;

    // Drms was already computed
    if (fabs(ediff) < energy_threshold_ && Drms_ < density_threshold_)
        return true;
    else
        return false;
}

void CUHF::form_initialF()
{
    Fa_->copy(H_);
    Fb_->copy(H_);

    if (debug_) {
        outfile->Printf( "Initial Fock alpha matrix:\n");
        Fa_->print("outfile");
        outfile->Printf( "Initial Fock beta matrix:\n");
        Fb_->print("outfile");
    }
}

void CUHF::form_F()
{
    // Form (rho_a + rho_b) / 2
    Dp_->copy(Da_);
    Dp_->add(Db_);
    Dp_->scale(-0.5); // This is a hack to get the eigenvectors in the
                      // order that I want
    if (debug_) {
      outfile->Printf( "Charge Density Matrix (SO Basis):\n");
      Dp_->print();
    }

    // Transfrom to an orthonormal basis, C_a is convenient
    Dp_->transform(S_);
    Dp_->transform(Ca_);
    if (debug_) {
      outfile->Printf( "Charge Density Matrix (Alpha Basis):\n");
      Dp_->print();
    }

    // Diagonalize the charge density and form the natural orbitals
    Dp_->diagonalize(Cno_temp_,No_);
    if (debug_) {
      outfile->Printf( "CUHF Natural Orbital Occupations:\n");
      No_->print();
    }
    Cno_->gemm(false, false, 1.0, Ca_, Cno_temp_, 0.0);

    // Now we form the contributions to the Fock matrix from
    // the charge and spin densities
    Fp_->copy(J_);
    Fp_->scale(2.0);
    Fp_->subtract(Ka_);
    Fp_->subtract(Kb_);
    Fp_->scale(0.5);

    Fm_->copy(Ka_);
    Fm_->subtract(Kb_);
    Fm_->scale(-0.5);

    // Transform the spin density contributions to the NO basis
    Fm_->transform(Cno_);

    // Zero the core-virtual contributions
    //
    //            [ Fm_cc Fm_co   0   ]
    // Fm_tilde = [ Fm_oc Fm_oo Fm_ov ]
    //            [   0   Fm_vo Fm_vv ]
    //
    for (int h = 0; h < nirrep_; ++h) {
      for (int i = 0; i < doccpi_[h]; ++i) {
        for (int j = doccpi_[h] + soccpi_[h]; j < nmopi_[h]; ++j) {
          Fm_->set(h, i, j, 0.0);
          Fm_->set(h, j, i, 0.0);
        }
      }
    }

    // Return to the SO basis
    Fm_->back_transform(Cno_);
    Fm_->transform(S_);

    // Build the modified alpha and beta Fock matrices
    Fa_->copy(H_);
    Fa_->add(Fp_);
    Fa_->add(Fm_);

    Fb_->copy(H_);
    Fb_->add(Fp_);
    Fb_->subtract(Fm_);

    if (debug_) {
        Fa_->print("outfile");
        Fb_->print("outfile");
    }
}

void CUHF::form_C()
{
    diagonalize_F(Fa_, Ca_, epsilon_a_);
    diagonalize_F(Fb_, Cb_, epsilon_b_);
    find_occupation();
    if (debug_) {
        Ca_->print("outfile");
        Cb_->print("outfile");
    }
}

void CUHF::form_D()
{
    for (int h = 0; h < nirrep_; ++h) {
        int nso = nsopi_[h];
        int nmo = nmopi_[h];
        int na = nalphapi_[h];
        int nb = nbetapi_[h];

        if (nso == 0 || nmo == 0) continue;

        double** Ca = Ca_->pointer(h);
        double** Cb = Cb_->pointer(h);
        double** Da = Da_->pointer(h);
        double** Db = Db_->pointer(h);

        if (na == 0)
            memset(static_cast<void*>(Da[0]), '\0', sizeof(double)*nso*nso);
        if (nb == 0)
            memset(static_cast<void*>(Db[0]), '\0', sizeof(double)*nso*nso);

        C_DGEMM('N','T',nso,nso,na,1.0,Ca[0],nmo,Ca[0],nmo,0.0,Da[0],nso);
        C_DGEMM('N','T',nso,nso,nb,1.0,Cb[0],nmo,Cb[0],nmo,0.0,Db[0],nso);

    }

    Dt_->copy(Da_);
    Dt_->add(Db_);

    if (debug_) {
        outfile->Printf( "in CUHF::form_D:\n");
        Da_->print();
        Db_->print();
    }
}

// TODO: Once Dt_ is refactored to D_ the only difference between this and RHF::compute_initial_E is a factor of 0.5
double CUHF::compute_initial_E()
{
    return nuclearrep_ + 0.5 * (Dt_->vector_dot(H_));
}

double CUHF::compute_E()
{
    double one_electron_E = Dt_->vector_dot(H_);
    double two_electron_E = 0.5 * (Da_->vector_dot(Fa_) + Db_->vector_dot(Fb_) - one_electron_E);

    energies_["Nuclear"] = nuclearrep_;
    energies_["One-Electron"] = one_electron_E;
    energies_["Two-Electron"] = two_electron_E;
    energies_["XC"] = 0.0;
    energies_["-D"] = 0.0;

    double DH  = Dt_->vector_dot(H_);
    double DFa = Da_->vector_dot(Fa_);
    double DFb = Db_->vector_dot(Fb_);
    double Eelec = 0.5 * (DH + DFa + DFb);
    // outfile->Printf( "electronic energy = %20.14f\n", Eelec);
    double Etotal = nuclearrep_ + Eelec;
    return Etotal;
}

void CUHF::compute_orbital_gradient(bool save_diis)
{
    SharedMatrix grad_a = form_FDSmSDF(Fa_, Da_);
    SharedMatrix grad_b = form_FDSmSDF(Fb_, Db_);

    // Store the RMS gradient for convergence checking
    Drms_ = 0.5 * (grad_a->rms() + grad_b->rms());

    if (save_diis){
        if (initialized_diis_manager_ == false) {
            diis_manager_ = std::shared_ptr<DIISManager>(new DIISManager(
                                                               max_diis_vectors_, "HF DIIS vector", DIISManager::LargestError,
                                                               DIISManager::OnDisk));
            diis_manager_->set_error_vector_size(2, DIISEntry::Matrix,
                                                 grad_a.get(), DIISEntry::Matrix, grad_b.get());
            diis_manager_->set_vector_size(2, DIISEntry::Matrix,
                                           Fa_.get(), DIISEntry::Matrix, Fb_.get());
            initialized_diis_manager_ = true;
        }

        diis_manager_->add_entry(4, grad_a.get(), grad_b.get(), Fa_.get(),
                                 Fb_.get());
    }
}

bool CUHF::diis()
{
    return diis_manager_->extrapolate(2, Fa_.get(), Fb_.get());
}

bool CUHF::stability_analysis()
{
    throw PSIEXCEPTION("CUHF stability analysis has not been implemented yet.  Sorry :(");
    return false;
}

}}
