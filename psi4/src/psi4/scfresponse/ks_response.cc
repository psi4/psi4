/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2026 The Psi4 Developers.
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

#include <cmath>
#include <string>

#include <psi4/psi4-dec.h>
#include <psi4/libpsi4util/process.h>
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

    // ensure running in c1 symmetry
    if ( reference_wavefunction_->nirrep() > 1 ) {
        throw PsiException("SCF response: only C1 symmetry is supported; add 'symmetry c1' to the molecule.",__FILE__,__LINE__);
    }

    // restricted vs unrestricted follows the reference orbitals, not the REFERENCE keyword
    // (a DFT calculation with the default REFERENCE RHF produces an RKS wavefunction)
    restricted_ = reference_wavefunction_->same_a_b_orbs();
    std::string reference = options_.get_str("REFERENCE");
    if ( reference == "ROHF" || reference == "ROKS" || reference == "CUHF" || (restricted_ && nalpha_ != nbeta_) ) {
        throw PsiException("SCF response: ROHF, ROKS, and CUHF references are not supported; use RHF/RKS or UHF/UKS.",__FILE__,__LINE__);
    }

    // get primary basis:
    std::shared_ptr<BasisSet> primary = reference_wavefunction_->get_basisset("ORBITAL");

    // determine the DFT functional and initialize the potential object. The potential is only
    // needed for the XC kernel; for Hartree-Fock Psi4 builds a polarized placeholder functional
    // even for RHF, which a restricted (RV) potential would reject.
    psi::scf::HF* scfwfn = (psi::scf::HF*)reference_wavefunction_.get();
    std::shared_ptr<SuperFunctional> functional = scfwfn->functional();

    // unsupported functional features: stop before any work (VBase::compute_Vx has no kernel for
    // them, and it throws only after turning a timer on, which breaks later calls in the process)
    if ( functional->is_meta() ) {
        throw PsiException("SCF response: meta-GGA functionals are not supported (no meta-GGA XC kernel).",__FILE__,__LINE__);
    }
    if ( functional->needs_vv10() ) {
        throw PsiException("SCF response: VV10 nonlocal correlation is not supported (no VV10 XC kernel).",__FILE__,__LINE__);
    }
    if ( options_.get_bool("WCOMBINE") ) {
        throw PsiException("SCF response: WCOMBINE = true is not supported.",__FILE__,__LINE__);
    }

    if ( functional->needs_xc() ) {
        potential_ = (std::shared_ptr<VBase>)VBase::build_V(primary,functional,options_,(restricted_ ? "RV" : "UV"));
        potential_->initialize();
    }

    is_hf_ = functional->name() == "HF" ? true : false;

    is_x_lrc_    = functional->is_x_lrc();
    is_x_hybrid_ = functional->is_x_hybrid();
    x_omega_     = functional->x_omega();
    if ( options_["DFT_OMEGA"].has_changed() ) {
        x_omega_ = options_.get_double("DFT_OMEGA");
    }
    x_alpha_     = functional->x_alpha();
    x_beta_      = functional->x_beta();
    needs_xc_    = functional->needs_xc();

    // JK object: reuse the one the SCF built (kept when SAVE_JK is on, which properties() sets),
    // so every SCF_TYPE works and the response uses the same integrals as the SCF. Otherwise
    // build one the same way the SCF does.
    if ( scfwfn->jk() ) {
        jk_ = scfwfn->jk();
    } else {
        std::shared_ptr<BasisSet> auxiliary = reference_wavefunction_->basisset_exists("DF_BASIS_SCF")
                                            ? reference_wavefunction_->get_basisset("DF_BASIS_SCF")
                                            : BasisSet::zero_ao_basis_set();
        // memory for jk in doubles (80% of what is available)
        size_t jk_doubles = (size_t)(0.8 * memory_ / sizeof(double));
        jk_ = JK::build_JK(primary, auxiliary, options_, is_x_lrc_, jk_doubles);
        jk_->set_memory(jk_doubles);
        jk_->set_cutoff(options_.get_double("INTS_TOLERANCE"));
        jk_->set_do_J(true);
        jk_->set_do_K(is_x_hybrid_);
        jk_->set_do_wK(is_x_lrc_);
        jk_->set_omega(x_omega_);
        jk_->initialize();
    }
    // the response products always need J
    jk_->set_do_J(true);

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

// Wavelength label used in variable names, as for the CC response properties:
// the wavelength in nm rounded to the nearest integer, or "INF NM" for a static field.
static std::string frequency_label(double omega) {
    if (std::fabs(omega) < 1.0e-14) return "INF NM";
    double nm = (pc_c * pc_h * 1e9) / (pc_hartree2J * std::fabs(omega));
    return std::to_string(std::lround(nm)) + "NM";
}

// Frequency line for the output file: omega in au and the wavelength in nm.
static void print_frequency(double omega) {
    if (std::fabs(omega) < 1.0e-14) {
        outfile->Printf("    omega = 0 (static field)\n\n");
    } else {
        double nm = (pc_c * pc_h * 1e9) / (pc_hartree2J * std::fabs(omega));
        outfile->Printf("    omega = %.8f au (%.2f nm)\n\n", omega, nm);
    }
}

/*- Process::environment.globals["SCF DIPOLE POLARIZABILITY @ xNM"] -*/
/*- Process::environment.arrays["SCF DIPOLE POLARIZABILITY TENSOR @ xNM"] -*/
void KSResponseSolver::store_polarizability(SharedMatrix alpha, double omega) {
    std::string at = " @ " + frequency_label(omega);
    double isotropic = alpha->trace() / 3.0;

    const char* xyz = "XYZ";
    outfile->Printf("\n    ==> Dipole Polarizability <==\n\n");
    print_frequency(omega);
    outfile->Printf("    %6s%20s%20s%20s\n", "", "X", "Y", "Z");
    for (int p = 0; p < 3; p++) {
        outfile->Printf("    %6c", xyz[p]);
        for (int q = 0; q < 3; q++) outfile->Printf("%20.12f", alpha->get(p, q));
        outfile->Printf("\n");
    }
    outfile->Printf("\n    Isotropic polarizability: %20.12f au\n", isotropic);

    Process::environment.globals["SCF DIPOLE POLARIZABILITY" + at] = isotropic;  // no-autodoc
    Process::environment.arrays["SCF DIPOLE POLARIZABILITY TENSOR" + at] = alpha->clone();  // no-autodoc
    reference_wavefunction_->set_scalar_variable("SCF DIPOLE POLARIZABILITY" + at, isotropic);  // no-autodoc
    reference_wavefunction_->set_array_variable("SCF DIPOLE POLARIZABILITY TENSOR" + at, alpha->clone());  // no-autodoc
}

/*- Process::environment.arrays["SCF DIPOLE HYPERPOLARIZABILITY TENSOR"] -*/
/*- Process::environment.arrays["SCF SHG HYPERPOLARIZABILITY TENSOR @ xNM"] -*/
/*- Process::environment.arrays["SCF OPTICAL RECTIFICATION HYPERPOLARIZABILITY TENSOR @ xNM"] -*/
/*- Process::environment.arrays["SCF POCKELS HYPERPOLARIZABILITY TENSOR @ xNM"] -*/
void KSResponseSolver::store_hyperpolarizability(SharedMatrix beta, const std::string& type, double omega) {
    std::string key;
    if (type == "STATIC") {
        key = "SCF DIPOLE HYPERPOLARIZABILITY TENSOR";
    } else if (type == "SHG") {
        key = "SCF SHG HYPERPOLARIZABILITY TENSOR @ " + frequency_label(omega);
    } else if (type == "OR") {
        key = "SCF OPTICAL RECTIFICATION HYPERPOLARIZABILITY TENSOR @ " + frequency_label(omega);
    } else if (type == "POCKELS") {
        key = "SCF POCKELS HYPERPOLARIZABILITY TENSOR @ " + frequency_label(omega);
    } else {
        throw PsiException("SCF response: unknown hyperpolarizability type " + type, __FILE__, __LINE__);
    }

    std::string title = type == "STATIC" ? "static, beta(0;0,0)"
                      : type == "SHG"    ? "second-harmonic generation, beta(-2w;w,w)"
                      : type == "OR"     ? "optical rectification, beta(0;w,-w)"
                                         : "Pockels effect, beta(-w;w,0)";
    const char* xyz = "XYZ";
    outfile->Printf("\n    ==> First Hyperpolarizability: %s <==\n\n", title.c_str());
    print_frequency(omega);
    outfile->Printf("    beta(i,j,k); i is the induced dipole, j and k the field directions\n\n");
    for (int p = 0; p < 3; p++) {
        outfile->Printf("    i = %c %20s%20s%20s\n", xyz[p], "k = X", "k = Y", "k = Z");
        for (int q = 0; q < 3; q++) {
            outfile->Printf("    j = %c ", xyz[q]);
            for (int r = 0; r < 3; r++) outfile->Printf("%20.12f", beta->get(p, 3 * q + r));
            outfile->Printf("\n");
        }
        outfile->Printf("\n");
    }
    Process::environment.arrays[key] = beta->clone();  // no-autodoc
    reference_wavefunction_->set_array_variable(key, beta->clone());
}

void KSResponseSolver::compute_properties(const std::string& property){

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
                throw PsiException("SCF response: unknown unit for OMEGA; use au, Hz, nm, or eV.", __FILE__, __LINE__);
            }
            omega.push_back(tmp_omega);
    
        }
    }

    // check the request before doing any work
    if ( property != "POLARIZABILITY" && !is_hf_ ) {
        throw PsiException("SCF response: hyperpolarizabilities are only available for Hartree-Fock; DFT would need "
                           "the third-order XC kernel, which is not implemented.", __FILE__, __LINE__);
    }
    if ( property == "HYPERPOLARIZABILITY" ) {
        for (auto my_omega: omega) {
            if (fabs(my_omega) > 1e-14) {
                throw PsiException("SCF response: SCF_DIPOLE_HYPERPOLARIZABILITY is the static hyperpolarizability, but OMEGA "
                                   "is nonzero. For frequency-dependent ones use SCF_SHG, SCF_OPTICAL_RECTIFICATION, or SCF_POCKELS.",
                                   __FILE__, __LINE__);
            }
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
    
    if ( property == "POLARIZABILITY" ) {
        for (auto my_omega: omega) {
            std::vector<std::vector<double>> amps = first_order_response(omega_a, omega_b, my_omega);
            compute_polarizability(amps[0], amps[1], my_omega);
        }
    }else if ( property == "HYPERPOLARIZABILITY" ) {
        for (auto my_omega: omega) {
            std::vector<std::vector<double>> amps = first_order_response(omega_a, omega_b, my_omega);
            compute_polarizability(amps[0], amps[1], my_omega);
            compute_hyperpolarizability(amps, amps, amps, "STATIC", my_omega);
        }
    }else if ( property == "SHG" ) {
        for (auto my_omega: omega) {
            std::vector<std::vector<double>> amps_pw = first_order_response(omega_a, omega_b, my_omega);
            compute_polarizability(amps_pw[0], amps_pw[1], my_omega);
            std::vector<std::vector<double>> amps_mtw = first_order_response(omega_a, omega_b, -2 * my_omega);
            compute_polarizability(amps_mtw[0], amps_mtw[1], -2 * my_omega);
            compute_hyperpolarizability(amps_mtw, amps_pw, amps_pw, "SHG", my_omega);
        }
    }else if ( property == "OR" ) {
        for (auto my_omega: omega) {
            std::vector<std::vector<double>> amps_0 = first_order_response(omega_a, omega_b, 0.0);
            compute_polarizability(amps_0[0], amps_0[1], 0.0);
            std::vector<std::vector<double>> amps_pw = first_order_response(omega_a, omega_b, my_omega);
            compute_polarizability(amps_pw[0], amps_pw[1], my_omega);
            std::vector<std::vector<double>> amps_mw = {amps_pw[1], amps_pw[0]};
            compute_hyperpolarizability(amps_0, amps_pw, amps_mw, "OR", my_omega);
        }
    }else if ( property == "POCKELS" ) {
        for (auto my_omega: omega) {
            std::vector<std::vector<double>> amps_0 = first_order_response(omega_a, omega_b, 0.0);
            compute_polarizability(amps_0[0], amps_0[1], 0.0);
            std::vector<std::vector<double>> amps_pw = first_order_response(omega_a, omega_b, my_omega);
            compute_polarizability(amps_pw[0], amps_pw[1], my_omega);
            std::vector<std::vector<double>> amps_mw = {amps_pw[1], amps_pw[0]};
            compute_hyperpolarizability(amps_mw, amps_pw, amps_0, "POCKELS", my_omega);
        }
    }else if ( property == "ALL" ) {
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
        throw PsiException("SCF response: unknown property type.",__FILE__,__LINE__);
    }
}

}} // End namespaces

