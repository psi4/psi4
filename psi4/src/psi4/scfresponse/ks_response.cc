/* 
 *  @BEGIN LICENSE
 * 
 *  Hilbert: a space for quantum chemistry plugins to Psi4 
 * 
 *  Copyright (c) 2020 by its authors (LICENSE).
 * 
 *  The copyrights for code used from other parties are included in
 *  the corresponding files.
 * 
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 * 
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 * 
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with this program.  If not, see http://www.gnu.org/licenses/.
 * 
 *  @END LICENSE
 */

#include <psi4/psi4-dec.h>
#include <psi4/physconst.h>
#include <psi4/liboptions/liboptions.h>
#include <psi4/libpsio/psio.hpp>

#include <psi4/libmints/molecule.h>
#include <psi4/libmints/wavefunction.h>
#include <psi4/libmints/mintshelper.h>
#include <psi4/libmints/matrix.h>
#include <psi4/libmints/factory.h>
#include <psi4/libmints/vector.h>
#include <psi4/libmints/basisset.h>
#include <psi4/lib3index/dftensor.h>
#include <psi4/libqt/qt.h>

// jk object
#include <psi4/libfock/jk.h>

// for dft
#include "psi4/libfock/v.h"
#include "psi4/libfunctional/superfunctional.h"
#include "psi4/libscf_solver/hf.h"

#include <psi4/psifiles.h>
#include <psi4/libtrans/integraltransform.h>

#include <algorithm>

// plugin headers
#include "uks_response.h"
#include "diis.h"

namespace psi { namespace scf_response{ 

KSResponseSolver::KSResponseSolver(std::shared_ptr<Wavefunction> reference_wavefunction, Options& options_):
    Wavefunction(options_) {
    reference_wavefunction_ = reference_wavefunction;
    common_init();
}

KSResponseSolver::~KSResponseSolver() {
}

void KSResponseSolver::common_init() {

    // from reference:

    shallow_copy(reference_wavefunction_);

    energy_   = reference_wavefunction_->energy();
    nalpha_   = reference_wavefunction_->nalpha();
    nbeta_    = reference_wavefunction_->nbeta();
    nalphapi_ = reference_wavefunction_->nalphapi();
    nbetapi_  = reference_wavefunction_->nbetapi();
    frzcpi_   = reference_wavefunction_->frzcpi();
    frzvpi_   = reference_wavefunction_->frzvpi();
    nmopi_    = reference_wavefunction_->nmopi();
    nirrep_   = reference_wavefunction_->nirrep();
    nso_      = reference_wavefunction_->nso();
    nmo_      = reference_wavefunction_->nmo();
    nsopi_    = reference_wavefunction_->nsopi();
    molecule_ = reference_wavefunction_->molecule();

    AO2SO_ = std::shared_ptr<Matrix>(reference_wavefunction_->aotoso());

    Ca_ = (std::shared_ptr<Matrix>)(new Matrix(reference_wavefunction_->Ca()));
    Cb_ = (std::shared_ptr<Matrix>)(new Matrix(reference_wavefunction_->Cb()));

    S_  = (std::shared_ptr<Matrix>)(new Matrix(reference_wavefunction_->S()));

    H_  = (std::shared_ptr<Matrix>)(new Matrix(reference_wavefunction_->H()));

    Fa_  = (std::shared_ptr<Matrix>)(new Matrix(reference_wavefunction_->Fa()));
    Fb_  = (std::shared_ptr<Matrix>)(new Matrix(reference_wavefunction_->Fb()));

    Da_  = (std::shared_ptr<Matrix>)(new Matrix(reference_wavefunction_->Da()));
    Db_  = (std::shared_ptr<Matrix>)(new Matrix(reference_wavefunction_->Db()));

    // Lagrangian matrix
    Lagrangian_ = std::shared_ptr<Matrix>(reference_wavefunction_->lagrangian());

    epsilon_a_ = std::make_shared<Vector>(nmopi_);
    epsilon_a_->copy(*reference_wavefunction_->epsilon_a().get());
    epsilon_b_ = std::make_shared<Vector>(nmopi_);
    epsilon_b_->copy(*reference_wavefunction_->epsilon_b().get());

    // other stuff:
    gradient_     =  reference_wavefunction_->matrix_factory()->create_shared_matrix("Total gradient", molecule_->natom(), 3);

    // memory is from process::environment
    memory_ = Process::environment.get_memory();

    // check SCF type
    if ( options_.get_str("SCF_TYPE") != "DF" && options_.get_str("SCF_TYPE") != "CD" && options_.get_str("SCF_TYPE") != "PK") {
        throw PsiException("invalid SCF_TYPE for scf_response",__FILE__,__LINE__);
    }

    // ensure running in c1 symmetry
    if ( reference_wavefunction_->nirrep() > 1 ) {
        throw PsiException("scf_response only works with c1 symmetry for now.",__FILE__,__LINE__);
    }

    // get primary basis:
    std::shared_ptr<BasisSet> primary = reference_wavefunction_->get_basisset("ORBITAL");

    // determine the DFT functional and initialize the potential object
    psi::scf::HF* scfwfn = (psi::scf::HF*)reference_wavefunction_.get();
    std::shared_ptr<SuperFunctional> functional = scfwfn->functional();
    potential_ = (std::shared_ptr<VBase>)VBase::build_V(primary,functional,options_,(options_.get_str("REFERENCE") == "RKS" ? "RV" : "UV"));

    // initialize potential object
    potential_->initialize();

    // apparently compute_Vx wants me to set the density (but this should be done in the uks/rks derived class)
    //potential_->set_D({Da_,Db_});

    if ( options_.get_str("SCF_TYPE") == "DF" ) {

        // get auxiliary basis:
        std::shared_ptr<BasisSet> auxiliary = reference_wavefunction_->get_basisset("DF_BASIS_SCF");

        jk_ = (std::shared_ptr<DiskDFJK>)(new DiskDFJK(primary,auxiliary,options_));

    }else if ( options_.get_str("SCF_TYPE") == "CD" ) {

        jk_ = (std::shared_ptr<CDJK>)(new CDJK(primary,options_,options_.get_double("CHOLESKY_TOLERANCE")));

    }else if ( options_.get_str("SCF_TYPE") == "PK" ) {

        jk_ = (std::shared_ptr<PKJK>)(new PKJK(primary,options_));
    }

    // memory for jk (say, 80% of what is available)
    jk_->set_memory(0.8 * memory_);

    // integral cutoff
    jk_->set_cutoff(options_.get_double("INTS_TOLERANCE"));

    is_hf_ = functional->name() == "HF" ? true : false;

    is_x_lrc_    = functional->is_x_lrc();
    is_x_hybrid_ = functional->is_x_hybrid();
    x_omega_     = functional->x_omega();
    if ( options_["DFT_OMEGA"].has_changed() ) {
        x_omega_ = options_.get_double("DFT_OMEGA");
    }
    x_alpha_     = functional->x_alpha();
    needs_xc_    = functional->needs_xc();

    jk_->set_do_J(true);
    jk_->set_do_K(is_x_hybrid_);
    jk_->set_do_wK(is_x_lrc_);
    jk_->set_omega(x_omega_);

    jk_->initialize();

    // dipole integrals
    std::shared_ptr<MintsHelper> mints (new MintsHelper(reference_wavefunction_));
    std::vector<std::shared_ptr<Matrix>> dipole = mints->so_dipole();

    mua_.clear();
    for (int i = 0; i < 3; i++) { 
        std::shared_ptr<Matrix> tmp = dipole[i]->clone();
        mua_.push_back(tmp);
        mua_[i]->transform(Ca_);
    }       
    mub_.clear();
    for (int i = 0; i < 3; i++) {
        std::shared_ptr<Matrix> tmp = dipole[i]->clone();
        mub_.push_back(tmp); 
        mub_[i]->transform(Cb_); 
    }
}

void KSResponseSolver::compute_properties(){

    int n_omega = 0;
    std::vector<double> omega;
    std::string units = "AU";
    
    // grab the field freqs from input -- a few units are converted to E_h
    int count = (int)options_["OMEGA"].size();
    if (count == 0) {  // Assume 0.0 E_h for field energy
        n_omega = 1;
        omega.push_back(0.0);
    } else if (count == 1) {  // Assume E_h for field energy and read value
        n_omega = 1;
        omega.push_back(options_["OMEGA"][0].to_double());
    } else if (count >= 2) {
        n_omega = count - 1;
        units = options_["OMEGA"][count - 1].to_string();
        for (int i = 0; i < count-1; i++) {
            double tmp_omega = options_["OMEGA"][i].to_double();
    
            if (units == "HZ" || units == "Hz" || units == "hz") {
                tmp_omega *= pc_h / pc_hartree2J;
            }else if (units == "AU" || units == "Au" || units == "au") { // do nothing
            }else if (units == "NM" || units == "nm") {
                tmp_omega = (pc_c * pc_h * 1e9) / (tmp_omega * pc_hartree2J);
            }else if (units == "EV" || units == "ev" || units == "eV") {
                tmp_omega /= pc_hartree2ev;
            }else {
                throw PsiException("Error in unit for input field frequencies, should be au, Hz, nm, or eV", __FILE__, __LINE__);
            }
            omega.push_back(tmp_omega);
    
        }
    }

    // omega = -mu
    std::vector<std::shared_ptr<Matrix>> omega_a;
    for (int i = 0; i < 3; i++) {
        std::shared_ptr<Matrix> tmp = mua_[i]->clone();
        omega_a.push_back(tmp);
        omega_a[i]->scale(-1);
    }
    std::vector<std::shared_ptr<Matrix>> omega_b;
    for (int i = 0; i < 3; i++) {
        std::shared_ptr<Matrix> tmp = mub_[i]->clone();
        omega_b.push_back(tmp);
        omega_b[i]->scale(-1);
    }
    
    if ( options_.get_str("PROPERTY") == "POLARIZABILITY" ) {
        for (auto my_omega: omega) {
            std::vector<std::vector<double>> amps = first_order_response(omega_a, omega_b, my_omega);
            compute_polarizability(amps[0], amps[1], my_omega);
        }
    }else if ( options_.get_str("PROPERTY") == "HYPERPOLARIZABILITY" ) {
        for (auto my_omega: omega) {
            if (fabs(my_omega) > 1e-14) {
                throw PsiException("For frequency-dependent hyperpolarizabilities, use PROPERTY = SHG, OR, or POCKELS", __FILE__, __LINE__);
            }
            std::vector<std::vector<double>> amps = first_order_response(omega_a, omega_b, my_omega);
            compute_polarizability(amps[0], amps[1], my_omega);
            compute_hyperpolarizability(amps, amps, amps, "STATIC", my_omega);
        }
    }else if ( options_.get_str("PROPERTY") == "SHG" ) {
        for (auto my_omega: omega) {
            std::vector<std::vector<double>> amps_pw = first_order_response(omega_a, omega_b, my_omega);
            compute_polarizability(amps_pw[0], amps_pw[1], my_omega);
            std::vector<std::vector<double>> amps_mtw = first_order_response(omega_a, omega_b, -2 * my_omega);
            compute_polarizability(amps_mtw[0], amps_mtw[1], -2 * my_omega);
            compute_hyperpolarizability(amps_mtw, amps_pw, amps_pw, "SHG", my_omega);
        }
    }else if ( options_.get_str("PROPERTY") == "OR" ) {
        for (auto my_omega: omega) {
            std::vector<std::vector<double>> amps_0 = first_order_response(omega_a, omega_b, 0.0);
            compute_polarizability(amps_0[0], amps_0[1], 0.0);
            std::vector<std::vector<double>> amps_pw = first_order_response(omega_a, omega_b, my_omega);
            compute_polarizability(amps_pw[0], amps_pw[1], my_omega);
            std::vector<std::vector<double>> amps_mw = {amps_pw[1], amps_pw[0]};
            compute_hyperpolarizability(amps_0, amps_pw, amps_mw, "OR", my_omega);
        }
    }else if ( options_.get_str("PROPERTY") == "POCKELS" ) {
        for (auto my_omega: omega) {
            std::vector<std::vector<double>> amps_0 = first_order_response(omega_a, omega_b, 0.0);
            compute_polarizability(amps_0[0], amps_0[1], 0.0);
            std::vector<std::vector<double>> amps_pw = first_order_response(omega_a, omega_b, my_omega);
            compute_polarizability(amps_pw[0], amps_pw[1], my_omega);
            std::vector<std::vector<double>> amps_mw = {amps_pw[1], amps_pw[0]};
            compute_hyperpolarizability(amps_mw, amps_pw, amps_0, "POCKELS", my_omega);
        }
    }else if ( options_.get_str("PROPERTY") == "ALL" ) {
        for (auto my_omega: omega) {
            if (fabs(my_omega) < 1e-14) {
                std::vector<std::vector<double>> amps = first_order_response(omega_a, omega_b, my_omega);
                compute_polarizability(amps[0], amps[1], my_omega);
                compute_hyperpolarizability(amps, amps, amps, "STATIC", my_omega);
            }else {
                std::vector<std::vector<double>> amps_0 = first_order_response(omega_a, omega_b, 0.0);
                compute_polarizability(amps_0[0], amps_0[1], 0.0);
                std::vector<std::vector<double>> amps_pw = first_order_response(omega_a, omega_b, my_omega);
                compute_polarizability(amps_pw[0], amps_pw[1], my_omega);
                std::vector<std::vector<double>> amps_mtw = first_order_response(omega_a, omega_b, -2 * my_omega);
                compute_polarizability(amps_mtw[0], amps_mtw[1], -2 * my_omega);
                std::vector<std::vector<double>> amps_mw = {amps_pw[1], amps_pw[0]};

                compute_hyperpolarizability(amps_0, amps_0, amps_0, "STATIC", 0.0);
                compute_hyperpolarizability(amps_0, amps_pw, amps_mw, "OR", my_omega);
                compute_hyperpolarizability(amps_mw, amps_pw, amps_0, "POCKELS", my_omega);
                compute_hyperpolarizability(amps_mtw, amps_pw, amps_pw, "SHG", my_omega);
            }
        }
    }else {
        throw PsiException("unsupported PROPERTY for scf_response",__FILE__,__LINE__);
    }
}

}} // End namespaces

