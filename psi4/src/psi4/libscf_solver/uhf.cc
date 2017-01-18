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
#include <tuple>

#include "psi4/libfunctional/superfunctional.h"
#include "psi4/libciomr/libciomr.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/libiwl/iwl.hpp"
#include "psi4/libqt/qt.h"

#include "psi4/libpsi4util/libpsi4util.h"
#include "psi4/libfock/v.h"
#include "psi4/libfock/jk.h"
#include "psi4/physconst.h"
#include "psi4/libtrans/integraltransform.h"
#include "psi4/libdpd/dpd.h"
#include "uhf.h"
#include "stability.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/factory.h"

namespace psi { namespace scf {

UHF::UHF(SharedWavefunction ref_wfn, std::shared_ptr<SuperFunctional> func)
    : HF(ref_wfn, func, Process::environment.options, PSIO::shared_object())
{
    common_init();
}

UHF::UHF(SharedWavefunction ref_wfn, std::shared_ptr<SuperFunctional> func,
         Options& options, std::shared_ptr<PSIO> psio)
    : HF(ref_wfn, func, options, psio)
{
    common_init();
}

UHF::~UHF()
{
}

void UHF::common_init()
{

    Drms_ = 0.0;
    // TODO: Move that to the base object
    step_scale_ = options_.get_double("FOLLOW_STEP_SCALE");
    step_increment_ = options_.get_double("FOLLOW_STEP_INCREMENT");

    Fa_     = SharedMatrix(factory_->create_matrix("F alpha"));
    Fb_     = SharedMatrix(factory_->create_matrix("F beta"));
    Da_     = SharedMatrix(factory_->create_matrix("SCF alpha density"));
    Db_     = SharedMatrix(factory_->create_matrix("SCF beta density"));
    Dt_     = SharedMatrix(factory_->create_matrix("D total"));
    Da_old_ = SharedMatrix(factory_->create_matrix("Old alpha SCF density"));
    Db_old_ = SharedMatrix(factory_->create_matrix("Old beta SCF density"));
    Dt_old_ = SharedMatrix(factory_->create_matrix("D total old"));
    Lagrangian_ = SharedMatrix(factory_->create_matrix("Lagrangian"));
    Ca_     = SharedMatrix(factory_->create_matrix("C alpha"));
    Cb_     = SharedMatrix(factory_->create_matrix("C beta"));
    Ga_     = SharedMatrix(factory_->create_matrix("G alpha"));
    Gb_     = SharedMatrix(factory_->create_matrix("G beta"));
    J_      = SharedMatrix(factory_->create_matrix("J total"));
    Ka_     = SharedMatrix(factory_->create_matrix("K alpha"));
    Kb_     = SharedMatrix(factory_->create_matrix("K beta"));
    wKa_    = SharedMatrix(factory_->create_matrix("wK alpha"));
    wKb_    = SharedMatrix(factory_->create_matrix("wK beta"));

    epsilon_a_ = SharedVector(factory_->create_vector());
    epsilon_b_ = SharedVector(factory_->create_vector());

    same_a_b_dens_ = false;
    same_a_b_orbs_ = false;
}

void UHF::finalize()
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
    Ga_.reset();
    Gb_.reset();

    compute_nos();

    HF::finalize();
}

void UHF::save_density_and_energy()
{
    Da_old_->copy(Da_);
    Db_old_->copy(Db_);
    Dt_old_->copy(Dt_);
    Eold_ = E_;
}
void UHF::form_V()
{
    // Push the C matrix on
    std::vector<SharedMatrix> & C = potential_->C();
    C.clear();
    C.push_back(Ca_subset("SO", "OCC"));
    C.push_back(Cb_subset("SO", "OCC"));

    // Run the potential object
    potential_->compute();

    // Pull the V matrices off
    const std::vector<SharedMatrix> & V = potential_->V();
    Va_ = V[0];
    Vb_ = V[1];
}
void UHF::form_G()
{
    if (functional_->needs_xc()) {
        timer_on("RKS: Form V");
        form_V();
        Ga_->copy(Va_);
        Gb_->copy(Vb_);
        timer_off("RKS: Form V");
    } else {
        Ga_->zero();
        Gb_->zero();
    }

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
    const std::vector<SharedMatrix> & wK = jk_->wK();
    J_->copy(J[0]);
    J_->add(J[1]);
    if (functional_->is_x_hybrid()) {
        Ka_ = K[0];
        Kb_ = K[1];
    }
    if (functional_->is_x_lrc()) {
        wKa_ = wK[0];
        wKb_ = wK[1];
    }
    Ga_->add(J_);
    Gb_->add(J_);

    double alpha = functional_->x_alpha();
    double beta = 1.0 - alpha;

    if (alpha != 0.0) {
        Ga_->axpy(-alpha, Ka_);
        Gb_->axpy(-alpha, Kb_);
    }
    else {
        Ka_->zero();
        Kb_->zero();
    }

    if (functional_->is_x_lrc()) {
        Ga_->axpy(-beta, wKa_);
        Gb_->axpy(-beta, wKb_);
    }
    else {
        wKa_->zero();
        wKb_->zero();
    }
}

void UHF::save_information()
{
}

bool UHF::test_convergency()
{
    double ediff = E_ - Eold_;

    // Drms was computed earlier
    if (fabs(ediff) < energy_threshold_ && Drms_ < density_threshold_)
        return true;
    else
        return false;
}

void UHF::form_initialF()
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

void UHF::form_F()
{
    Fa_->copy(H_);
    Fa_->add(Ga_);

    Fb_->copy(H_);
    Fb_->add(Gb_);

    if (debug_) {
        Fa_->print("outfile");
        Fb_->print("outfile");
    }
}

void UHF::form_C()
{
    diagonalize_F(Fa_, Ca_, epsilon_a_);
    diagonalize_F(Fb_, Cb_, epsilon_b_);
    if (options_.get_bool("GUESS_MIX") && (iteration_ == 0)){
        if (Ca_->nirrep() == 1){
            outfile->Printf("  Mixing alpha HOMO/LUMO orbitals (%d,%d)\n\n",nalpha_,nalpha_ + 1);
            Ca_->rotate_columns(0,nalpha_ - 1,nalpha_, pc_pi * 0.25);
            Cb_->rotate_columns(0,nbeta_ - 1,nbeta_,-pc_pi * 0.25);
        }else{
            throw InputException("Warning: cannot mix alpha HOMO/LUMO orbitals. Run in C1 symmetry.", "to 'symmetry c1'", __FILE__, __LINE__);
        }
    }
    find_occupation();
    if (debug_) {
        Ca_->print("outfile");
        Cb_->print("outfile");
    }
}

void UHF::form_D()
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
            ::memset(static_cast<void*>(Da[0]), '\0', sizeof(double)*nso*nso);
        if (nb == 0)
            ::memset(static_cast<void*>(Db[0]), '\0', sizeof(double)*nso*nso);

        C_DGEMM('N','T',nso,nso,na,1.0,Ca[0],nmo,Ca[0],nmo,0.0,Da[0],nso);
        C_DGEMM('N','T',nso,nso,nb,1.0,Cb[0],nmo,Cb[0],nmo,0.0,Db[0],nso);

    }

    Dt_->copy(Da_);
    Dt_->add(Db_);

    if (debug_) {
        outfile->Printf( "in UHF::form_D:\n");
        Da_->print();
        Db_->print();
    }
}

// TODO: Once Dt_ is refactored to D_ the only difference between this and RHF::compute_initial_E is a factor of 0.5
double UHF::compute_initial_E()
{
    Dt_->copy(Da_);
    Dt_->add(Db_);
    return nuclearrep_ + 0.5 * (Dt_->vector_dot(H_));
}

double UHF::compute_E()
{
    // E_DFT = 2.0 D*H + D*J - \alpha D*K + E_xc
    double one_electron_E = Da_->vector_dot(H_);
    one_electron_E += Db_->vector_dot(H_);
    double coulomb_E = Da_->vector_dot(J_);
    coulomb_E += Db_->vector_dot(J_);

    double XC_E = 0.0;
    if (functional_->needs_xc()) {
        XC_E = potential_->quadrature_values()["FUNCTIONAL"];
    }

    double exchange_E = 0.0;
    double alpha = functional_->x_alpha();
    double beta = 1.0 - alpha;
    if (functional_->is_x_hybrid()) {
        exchange_E -= alpha * Da_->vector_dot(Ka_);
        exchange_E -= alpha * Db_->vector_dot(Kb_);
    }
    if (functional_->is_x_lrc()) {
        exchange_E -= beta * Da_->vector_dot(wKa_);
        exchange_E -= beta * Db_->vector_dot(wKb_);
    }

    energies_["Nuclear"] = nuclearrep_;
    energies_["One-Electron"] = one_electron_E;
    energies_["Two-Electron"] = 0.5 * (coulomb_E + exchange_E);
    energies_["XC"] = XC_E;
    energies_["-D"] = variables_["-D Energy"];
    double dashD_E = energies_["-D"];

    double Etotal = 0.0;
    Etotal += nuclearrep_;
    Etotal += one_electron_E;
    Etotal += 0.5 * coulomb_E;
    Etotal += 0.5 * exchange_E;
    Etotal += XC_E;
    Etotal += dashD_E;

    return Etotal;
}
void UHF::Hx(SharedMatrix x_a, SharedMatrix IFock_a, SharedMatrix Cocc_a,
             SharedMatrix Cvir_a, SharedMatrix ret_a,
             SharedMatrix x_b, SharedMatrix IFock_b, SharedMatrix Cocc_b,
             SharedMatrix Cvir_b, SharedMatrix ret_b)
{
    if (functional_->needs_xc()){
        throw PSIEXCEPTION("SCF: Cannot yet compute DFT Hessian-vector products.\n");
    }

    // => Effective one electron part <= //
    Dimension virpi_a = Cvir_a->colspi();
    Dimension virpi_b = Cvir_b->colspi();
    for (size_t h=0; h<nirrep_; h++){
        // Alpha
        if (nalphapi_[h] && virpi_a[h]){
            double** IFp = IFock_a->pointer(h);
            double** retp = ret_a->pointer(h);
            double** xp = x_a->pointer(h);

            // ret_ia = F_ij X_ja
            C_DGEMM('N','N',nalphapi_[h],virpi_a[h],nalphapi_[h],1.0,
                    IFp[0],nmopi_[h],
                    xp[0],virpi_a[h],0.0,retp[0],virpi_a[h]);

            // ret_ia -= X_ib F_ba
            C_DGEMM('N','N',nalphapi_[h],virpi_a[h],virpi_a[h],-1.0,
                    xp[0],virpi_a[h],
                    (IFp[nalphapi_[h]]+nalphapi_[h]),nmopi_[h],1.0,retp[0],virpi_a[h]);

        }
        // Beta
        if (nbetapi_[h] && virpi_b[h]){
            double** IFp = IFock_b->pointer(h);
            double** retp = ret_b->pointer(h);
            double** xp = x_b->pointer(h);

            // ret_ia = F_ij X_ja
            C_DGEMM('N','N',nbetapi_[h],virpi_b[h],nbetapi_[h],1.0,
                    IFp[0],nmopi_[h],
                    xp[0],virpi_b[h],0.0,retp[0],virpi_b[h]);

            // ret_ia -= X_ib F_ba
            C_DGEMM('N','N',nbetapi_[h],virpi_b[h],virpi_b[h],-1.0,
                    xp[0],virpi_b[h],
                    (IFp[nbetapi_[h]]+nbetapi_[h]),nmopi_[h],1.0,retp[0],virpi_b[h]);

        }
    }

    // => Two electron part <= //
    std::vector<SharedMatrix>& Cl = jk_->C_left();
    std::vector<SharedMatrix>& Cr = jk_->C_right();
    Cl.clear();
    Cr.clear();

    Cl.push_back(Cocc_a);
    Cl.push_back(Cocc_b);

    SharedMatrix R_a = Matrix::doublet(Cvir_a, x_a, false, true);
    SharedMatrix R_b = Matrix::doublet(Cvir_b, x_b, false, true);
    R_a->scale(-1.0);
    R_b->scale(-1.0);
    Cr.push_back(R_a);
    Cr.push_back(R_b);

    jk_->compute();

    // Just in case someone only clears out Cleft and gets very strange errors
    Cl.clear();
    Cr.clear();

    const std::vector<SharedMatrix>& J = jk_->J();
    const std::vector<SharedMatrix>& K = jk_->K();

    // D_nm = Cocc_ni x_ia Cvir_ma
    // Cocc_ni (4 * J[D]_nm - K[D]_nm - K[D]_mn) C_vir_ma
    J[0]->add(J[1]);
    J[0]->scale(2.0);
    J[1]->copy(J[0]);

    J[0]->subtract(K[0]);
    J[0]->subtract(K[0]->transpose());
    R_a->gemm(false, false, 1.0, J[0], Cocc_a, 0.0);
    ret_a->gemm(true, false, 1.0, R_a, Cvir_a, 1.0);
    ret_a->scale(-4.0);

    J[1]->subtract(K[1]);
    J[1]->subtract(K[1]->transpose());
    R_b->gemm(false, false, 1.0, J[1], Cocc_b, 0.0);
    ret_b->gemm(true, false, 1.0, R_b, Cvir_b, 1.0);
    ret_b->scale(-4.0);

    // Cleaup
    R_a.reset();
    R_b.reset();

}
void UHF::damp_update()
{
  Da_->scale(1.0 - damping_percentage_);
  Da_->axpy(damping_percentage_, Da_old_);
  Db_->scale(1.0 - damping_percentage_);
  Db_->axpy(damping_percentage_, Db_old_);
  Dt_->copy(Da_);
  Dt_->add(Db_);
}
int UHF::soscf_update(void)
{
    time_t start, stop;
    start = time(NULL);

    // => Build gradient and preconditioner <= //

    // Grab occ and vir orbitals
    SharedMatrix Cocc_a = Ca_subset("SO", "OCC");
    SharedMatrix Cvir_a = Ca_subset("SO", "VIR");
    Dimension virpi_a = Cvir_a->colspi();

    SharedMatrix Cocc_b = Cb_subset("SO", "OCC");
    SharedMatrix Cvir_b = Cb_subset("SO", "VIR");
    Dimension virpi_b = Cvir_b->colspi();

    // MO Fock Matrix (Inactive Fock in Helgaker's language)
    SharedMatrix IFock_a = Matrix::triplet(Ca_, Fa_, Ca_, true, false, false);
    SharedMatrix Gradient_a = SharedMatrix(new Matrix("Alpha Gradient",
                                           nirrep_, nalphapi_, virpi_a));
    SharedMatrix Precon_a = SharedMatrix(new Matrix("Alpha Precon",
                                         nirrep_, nalphapi_, virpi_a));

    SharedMatrix IFock_b = Matrix::triplet(Cb_, Fb_, Cb_, true, false, false);
    SharedMatrix Gradient_b = SharedMatrix(new Matrix("Beta Gradient",
                                           nirrep_, nbetapi_, virpi_b));
    SharedMatrix Precon_b = SharedMatrix(new Matrix("Beta Precon",
                                         nirrep_, nbetapi_, virpi_b));

    for (size_t h=0; h<nirrep_; h++){

        if (nalphapi_[h] && virpi_a[h]){
            double* gp = Gradient_a->pointer(h)[0];
            double* denomp = Precon_a->pointer(h)[0];
            double** fp = IFock_a->pointer(h);

            for (size_t i=0, target=0; i<nalphapi_[h]; i++){
                for (size_t a=nalphapi_[h]; a < nmopi_[h]; a++){
                    gp[target] = -4.0 * fp[i][a];
                    denomp[target++] = -4.0 * (fp[i][i] - fp[a][a]);
                }
            }
        }
        if (nbetapi_[h] && virpi_b[h]){
            double* gp = Gradient_b->pointer(h)[0];
            double* denomp = Precon_b->pointer(h)[0];
            double** fp = IFock_b->pointer(h);

            for (size_t i=0, target=0; i<nbetapi_[h]; i++){
                for (size_t a=nbetapi_[h]; a < nmopi_[h]; a++){
                    gp[target] = -4.0 * fp[i][a];
                    denomp[target++] = -4.0 * (fp[i][i] - fp[a][a]);
                }
            }
        }
    }

    // Make sure the MO gradient is reasonably small
    if ((Gradient_a->absmax() > 0.3) || (Gradient_b->absmax() > 0.3)){
        if (print_ > 1){
            outfile->Printf("    Gradient element too large for SOSCF, using DIIS.\n");
        }
        return 0;
    }

    if (soscf_print_){
        outfile->Printf("\n");
        outfile->Printf("    ==> SOUHF Iterations <==\n");
        outfile->Printf("    Maxiter     = %11d\n", soscf_max_iter_);
        outfile->Printf("    Miniter     = %11d\n", soscf_min_iter_);
        outfile->Printf("    Convergence = %11.3E\n", soscf_conv_);
        outfile->Printf("    ---------------------------------------\n");
        outfile->Printf("    %-4s   %11s     %10s\n", "Iter", "Residual RMS", "Time [s]");
        outfile->Printf("    ---------------------------------------\n");
    }

    // => Initial CG guess <= //
    SharedMatrix x_a = Gradient_a->clone();
    x_a->apply_denominator(Precon_a);
    SharedMatrix x_b = Gradient_b->clone();
    x_b->apply_denominator(Precon_b);

    // Calc hessian vector product, find residual and conditioned residual
    SharedMatrix r_a = Gradient_a->clone();
    SharedMatrix Ap_a = SharedMatrix(new Matrix("Ap_a", nirrep_, nalphapi_, virpi_a));
    SharedMatrix r_b = Gradient_b->clone();
    SharedMatrix Ap_b = SharedMatrix(new Matrix("Ap_b", nirrep_, nbetapi_, virpi_b));

    Hx(x_a, IFock_a, Cocc_a, Cvir_a, Ap_a,
       x_b, IFock_b, Cocc_b, Cvir_b, Ap_b);
    r_a->subtract(Ap_a);
    r_b->subtract(Ap_b);

    // Print iteration 0 timings and rms
    double rconv = r_a->sum_of_squares() + r_b->sum_of_squares();
    double grad_rms = Gradient_a->sum_of_squares() + Gradient_b->sum_of_squares();
    if (grad_rms < 1.e-14){
        grad_rms = 1.e-14; // Prevent rel denom from being too small
    }
    double rms = 0.5 * sqrt(rconv / grad_rms);
    stop = time(NULL);
    if (soscf_print_){
        outfile->Printf("    %-5s %11.3E %10ld\n", "Guess", rms, stop-start);
    }

    // Build new p and z vectors
    SharedMatrix z_a = r_a->clone();
    z_a->apply_denominator(Precon_a);
    SharedMatrix p_a = z_a->clone();

    SharedMatrix z_b = r_b->clone();
    z_b->apply_denominator(Precon_b);
    SharedMatrix p_b = z_b->clone();

    // => CG iterations <= //
    int fock_builds = 1;
    for (int cg_iter=1; cg_iter<soscf_max_iter_; cg_iter++) {

        // Calc hessian vector product
        Hx(p_a, IFock_a, Cocc_a, Cvir_a, Ap_a,
           p_b, IFock_b, Cocc_b, Cvir_b, Ap_b);
        fock_builds += 1;

        // Find factors and scale
        double rzpre = r_a->vector_dot(z_a) + r_b->vector_dot(z_b);
        double alpha = rzpre / (p_a->vector_dot(Ap_a) + p_b->vector_dot(Ap_b));
        if (std::isnan(alpha)){
            outfile->Printf("UHF::SOSCF Warning CG alpha is zero/nan. Stopping CG.\n");
            alpha = 0.0;
        }

        x_a->axpy(alpha, p_a);
        r_a->axpy(-alpha, Ap_a);

        x_b->axpy(alpha, p_b);
        r_b->axpy(-alpha, Ap_b);

        // Get residual
        double rconv = r_a->sum_of_squares() + r_b->sum_of_squares();
        double rms = 0.5 * sqrt(rconv / grad_rms);
        stop = time(NULL);
        if (soscf_print_){
            outfile->Printf("    %-5d %11.3E %10ld\n", cg_iter, rms, stop-start);
        }

        // Check convergence
        if (((rms < soscf_conv_) && (cg_iter >= soscf_min_iter_)) || (alpha==0.0)) {
            break;
        }

        // Update p and z
        z_a->copy(r_a);
        z_a->apply_denominator(Precon_a);

        z_b->copy(r_b);
        z_b->apply_denominator(Precon_b);

        double beta = (r_a->vector_dot(z_a) + r_b->vector_dot(z_b)) / rzpre;

        p_a->scale(beta);
        p_a->add(z_a);

        p_b->scale(beta);
        p_b->add(z_b);

    }
    if (soscf_print_){
        outfile->Printf("    ---------------------------------------\n");
        outfile->Printf("\n");
    }

    // => Rotate orbitals <= //
    rotate_orbitals(Ca_, x_a);
    rotate_orbitals(Cb_, x_b);

    // => Cleanup <= //
    Cocc_a.reset();     Cocc_b.reset();
    Cvir_a.reset();     Cvir_b.reset();
    IFock_a.reset();    IFock_b.reset();
    Precon_a.reset();   Precon_b.reset();
    Gradient_a.reset(); Gradient_b.reset();
    Ap_a.reset();       Ap_b.reset();
    z_a.reset();        z_b.reset();
    r_a.reset();        r_b.reset();
    p_a.reset();        p_b.reset();

    return fock_builds;
}

void UHF::compute_orbital_gradient(bool save_fock)
{
    SharedMatrix gradient_a = form_FDSmSDF(Fa_, Da_);
    SharedMatrix gradient_b = form_FDSmSDF(Fb_, Db_);
    Drms_ = 0.5*(gradient_a->rms() + gradient_b->rms());

    if(save_fock){
        if (initialized_diis_manager_ == false) {
            diis_manager_ = std::shared_ptr<DIISManager>(new DIISManager(max_diis_vectors_, "HF DIIS vector", DIISManager::LargestError, DIISManager::OnDisk));
            diis_manager_->set_error_vector_size(2,
                                                 DIISEntry::Matrix, gradient_a.get(),
                                                 DIISEntry::Matrix, gradient_b.get());
            diis_manager_->set_vector_size(2,
                                           DIISEntry::Matrix, Fa_.get(),
                                           DIISEntry::Matrix, Fb_.get());
            initialized_diis_manager_ = true;
        }

        diis_manager_->add_entry(4, gradient_a.get(), gradient_b.get(), Fa_.get(), Fb_.get());
    }
}

bool UHF::diis()
{
    return diis_manager_->extrapolate(2, Fa_.get(), Fb_.get());
}

bool UHF::stability_analysis()
{
    std::shared_ptr<UStab> stab = std::shared_ptr<UStab>(new UStab(shared_from_this(), options_));
    stab->compute_energy();
    SharedMatrix eval_sym = stab->analyze();
    outfile->Printf( "    Lowest UHF->UHF stability eigenvalues: \n");
    std::vector < std::pair < double,int > >  eval_print;
    for (int h = 0; h < eval_sym->nirrep(); ++h) {
        for (int i = 0; i < eval_sym->rowdim(h); ++i) {
            eval_print.push_back(std::make_pair(eval_sym->get(h,i,0),h));
        }
    }
    print_stability_analysis(eval_print);

    // And now, export the eigenvalues to a PSI4 array, mainly for testing purposes

    Process::environment.arrays["SCF STABILITY EIGENVALUES"] = eval_sym;
    if (stab->is_unstable() && options_.get_str("STABILITY_ANALYSIS") == "FOLLOW") {
        if (attempt_number_ == 1 ) {
            stab_val = stab->get_eigval();
        } else if (stab_val - stab->get_eigval() < 1e-4) {
            // We probably fell on the same minimum, increase step_scale_
            outfile->Printf("    Negative eigenvalue similar to previous one, wavefunction\n");
            outfile->Printf("    likely to be in the same minimum.\n");
            step_scale_ += step_increment_;
            outfile->Printf("    Modifying FOLLOW_STEP_SCALE to %f.\n", step_scale_);
        } else {
            stab_val = stab->get_eigval();
        }
       //     outfile->Printf( "OLD ORBS");
       //     Ca_->print();
        stab->rotate_orbs(step_scale_);
       //     outfile->Printf( "NEW ORBS");
       //     Ca_->print();

       // Ask politely SCF control for a new set of iterations
       return true;
    } else {
        outfile->Printf("    Stability analysis over.\n");
        // We are done, no more iterations
        return false;
    }

}

void UHF::compute_nos()
{
    // Compute UHF NOs and NOONs [J. Chem. Phys. 88, 4926 (1988)] -- TDC, 8/15

    // Build S^1/2
    SharedMatrix SHalf = S_->clone();
    SHalf->power(0.5);

    // Diagonalize S^1/2 Dt S^1/2
    SharedMatrix SDS = factory_->create_shared_matrix("S^1/2 Dt S^1/2");
    SDS->copy(Da_);
    SDS->add(Db_);
    SDS->transform(SHalf);

    SharedMatrix UHF_NOs = factory_->create_shared_matrix("UHF NOs");
    SharedVector UHF_NOONs(factory_->create_vector());
    SDS->diagonalize(UHF_NOs, UHF_NOONs, descending);

    // Print the NOONs -- code ripped off from OEProp::compute_no_occupations()
    int max_num;
    if(options_.get_str("PRINT_NOONS") == "ALL") max_num = nmo_;
    else max_num = to_integer(options_.get_str("PRINT_NOONS"));

    std::vector<std::tuple<double, int, int> > metric;
    for (int h = 0; h < UHF_NOONs->nirrep(); h++)
      for (int i = 0; i < UHF_NOONs->dimpi()[h]; i++)
        metric.push_back(std::tuple<double,int,int>(UHF_NOONs->get(h,i), h ,i));

    std::sort(metric.begin(), metric.end(), std::greater<std::tuple<double,int,int> >());
    int offset = nalpha_;
    int start_occ = offset - max_num;
    start_occ = (start_occ < 0 ? 0 : start_occ);
    int stop_vir = offset + max_num + 1;
    stop_vir = (int)((size_t)stop_vir >= metric.size() ? metric.size() : stop_vir);
    char** labels = basisset_->molecule()->irrep_labels();
    outfile->Printf( "\n  UHF NO Occupations:\n");
    for (int index = start_occ; index < stop_vir; index++) {
      if (index < offset) {
        outfile->Printf( "  HONO-%-2d: %4d%3s %9.7f\n", offset- index - 1,
        std::get<2>(metric[index])+1,labels[std::get<1>(metric[index])],
        std::get<0>(metric[index]));
      }
      else {
        outfile->Printf( "  LUNO+%-2d: %4d%3s %9.7f\n", index - offset,
        std::get<2>(metric[index])+1,labels[std::get<1>(metric[index])],
        std::get<0>(metric[index]));
      }
    }
    outfile->Printf( "\n");

    if(options_.get_bool("SAVE_UHF_NOS")){
        // Save the NOs to Ca and Cb. The resulting orbitals will be restricted.

        outfile->Printf( "  Saving the UHF Natural Orbitals.\n");

        SharedMatrix SHalf_inv = S_->clone();
        SHalf_inv->power(-0.5);
        Ca_->gemm(false, false, 1.0,SHalf_inv,UHF_NOs, 0.0);

        double actv_threshold = 0.02;

        // Transform the average Fock matrix to the NO basis
        SharedMatrix F_UHF_NOs = factory_->create_shared_matrix("Fock Matrix");
        F_UHF_NOs->copy(Fa_);
        F_UHF_NOs->add(Fb_);
        F_UHF_NOs->transform(Ca_);

        // Sort orbitals according to type (core,active,virtual) and energy
        std::vector<std::tuple<int, double, int, int> > sorted_nos;
        for (int h = 0; h < UHF_NOONs->nirrep(); h++){
            for (int i = 0; i < UHF_NOONs->dimpi()[h]; i++){
                double noon = UHF_NOONs->get(h,i);
                int type = 0; // core      NO >= 1.98
                if (noon < actv_threshold){
                    type = 2; // virtual   NO < 0.02
                }else if (noon < 2.0 - actv_threshold){
                    type = 1; // active    0.02 <= NO < 1.98
                }
                double epsilon = F_UHF_NOs->get(h,i,i);
                sorted_nos.push_back(std::tuple<int,double,int,int>(type,epsilon,h,i));
            }
        }
        std::sort(sorted_nos.begin(), sorted_nos.end());

        // Build the final set of UHF NOs
        std::vector<int> irrep_count(nirrep_,0);

        for (size_t i = 0; i < sorted_nos.size(); i++){
            int h = std::get<2>(sorted_nos[i]);
            int Ca_p = std::get<3>(sorted_nos[i]);
            int Cb_p = irrep_count[h];
            for (int mu = 0; mu < Ca_->colspi(h); mu++){
                double value = Ca_->get(h,mu,Ca_p);
                Cb_->set(h,mu,Cb_p,value);
            }
            irrep_count[h] += 1;
        }

        // Copy sorted orbitals to Ca
        Ca_->copy(Cb_);

        // Suggest an active space
        Dimension corepi(nirrep_);
        Dimension actvpi(nirrep_);
        for (size_t i = 0; i < sorted_nos.size(); i++){
            int type = std::get<0>(sorted_nos[i]);
            int h = std::get<2>(sorted_nos[i]);
            if (type == 0) corepi[h] += 1;
            if (type == 1) actvpi[h] += 1;
        }
        outfile->Printf("\n  Active Space from UHF-NOs (NO threshold = %.4f):\n\n",actv_threshold);

        outfile->Printf("    restricted_docc = [");
        for (int h = 0; h < nirrep_; h++){
            outfile->Printf("%s%d",h ? "," : "",corepi[h]);
        }
        outfile->Printf("]\n");

        outfile->Printf("    active = [");
        for (int h = 0; h < nirrep_; h++){
            outfile->Printf("%s%d",h ? "," : "",actvpi[h]);
        }
        outfile->Printf("]\n");
    }
}

}}
