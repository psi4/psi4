/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2024 The Psi4 Developers.
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

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <iostream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "psi4/psifiles.h"
#include "psi4/physconst.h"

#include "psi4/libciomr/libciomr.h"
#include "psi4/libpsio/psio.h"
#include "psi4/libiwl/iwl.hpp"
#include "psi4/libqt/qt.h"
#include "psi4/libfock/jk.h"
#include "psi4/libfock/v.h"
#include "psi4/libfunctional/superfunctional.h"

#include "psi4/libpsi4util/libpsi4util.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/mintshelper.h"
#include "psi4/libmints/petitelist.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/extern.h"
#include "psi4/libmints/factory.h"
#include "psi4/libmints/pointgrp.h"
#include "psi4/libmints/oeprop.h"
#include "psi4/libmints/orthog.h"
#include "psi4/libmints/integral.h"
#include "psi4/libmints/quadrupole.h"
#include "psi4/libmints/sobasis.h"

#include "hf.h"

#include "psi4/psi4-dec.h"

#ifdef USING_BrianQC

#include <use_brian_wrapper.h>
#include <brian_macros.h>
#include <brian_scf.h>

extern void checkBrian();
extern BrianCookie brianCookie;
extern bool brianEnable;

#endif

namespace psi {
namespace scf {

HF::HF(SharedWavefunction ref_wfn, std::shared_ptr<SuperFunctional> func, Options& options, std::shared_ptr<PSIO> psio)
    : Wavefunction(options), functional_(func) {
    shallow_copy(ref_wfn);
    psio_ = psio;
    common_init();
}

HF::~HF() {}

void HF::common_init() {
    attempt_number_ = 1;
    reset_occ_ = false;
    sad_ = false;
    module_ = "scf";
    frac_performed_ = false;

    // This quantity is needed fairly soon
    nirrep_ = factory_->nirrep();

    integral_threshold_ = options_.get_double("INTS_TOLERANCE");

    scf_type_ = options_.get_str("SCF_TYPE");

    H_ = factory_->create_matrix("One-electron Hamiltonian");
    X_ = factory_->create_matrix("X");

    // nmo_ and nmopi_ not determined at present.
    nso_ = 0;
    const Dimension& dimpi = factory_->colspi();
    for (int h = 0; h < factory_->nirrep(); h++) {
        nsopi_[h] = dimpi[h];
        nso_ += nsopi_[h];
    }

    energies_["Total Energy"] = 0.0;

    // Read in DOCC and SOCC from memory
    input_docc_ = false;
    Dimension docc(nirrep_);
    if (options_["DOCC"].has_changed()) {
        input_docc_ = true;
        // Map the symmetry of the input DOCC, to account for displacements
        auto ps = options_.get_str("PARENT_SYMMETRY");
        if (ps != "") {
            auto old_pg = std::make_shared<PointGroup>(ps);
            // This is one of a series of displacements;  check the dimension against the parent
            // point group
            size_t full_nirreps = old_pg->char_table().nirrep();
            if (options_["DOCC"].size() != full_nirreps)
                throw PSIEXCEPTION("Input DOCC array has the wrong dimensions");
            Dimension temp_docc(full_nirreps);
            for (int h = 0; h < full_nirreps; ++h) {
                temp_docc[h] = options_["DOCC"][h].to_integer();
            }
            docc = map_irreps(temp_docc);
        } else {
            // This is a normal calculation; check the dimension against the current point group
            // then read
            if (options_["DOCC"].size() != nirrep_) throw PSIEXCEPTION("Input DOCC array has the wrong dimensions");
            for (int h = 0; h < nirrep_; ++h) {
                docc[h] = options_["DOCC"][h].to_integer();
            }
        }
    }  // else take the reference wavefunctions doccpi

    input_socc_ = false;
    Dimension socc(nirrep_);
    if (options_["SOCC"].has_changed()) {
        input_socc_ = true;
        // Map the symmetry of the input SOCC, to account for displacements
        auto ps = options_.get_str("PARENT_SYMMETRY");
        if (ps != "") {
            auto old_pg = std::make_shared<PointGroup>(ps);
            // This is one of a series of displacements;  check the dimension against the parent
            // point group
            size_t full_nirreps = old_pg->char_table().nirrep();
            if (options_["SOCC"].size() != full_nirreps)
                throw PSIEXCEPTION("Input SOCC array has the wrong dimensions");
            Dimension temp_socc(full_nirreps);
            for (int h = 0; h < full_nirreps; ++h) {
                temp_socc[h] = options_["SOCC"][h].to_integer();
            }
            socc = map_irreps(temp_socc);
        } else {
            // This is a normal calculation; check the dimension against the current point group
            // then read
            if (options_["SOCC"].size() != nirrep_) throw PSIEXCEPTION("Input SOCC array has the wrong dimensions");
            for (int h = 0; h < nirrep_; ++h) {
                socc[h] = options_["SOCC"][h].to_integer();
            }
        }
    }  // else take the reference wavefunctions soccpi

    if (input_socc_ || input_docc_) {
        nalphapi_ = docc + socc;
        nbetapi_ = docc;
        int alphacount = nalphapi_.sum();
        int betacount = nbetapi_.sum();
        if (alphacount != nalpha_) {
            std::ostringstream oss;
            oss << "Got " << alphacount << " alpha electrons, expected " << nalpha_ << ".\n";
            oss << "DOCC and SOCC must specify the occupation of all electrons or none.";
            throw PSIEXCEPTION(oss.str());
        }
        if (betacount != nbeta_) {
            std::ostringstream oss;
            oss << "Got " << betacount << " beta electrons, expected " << nbeta_ << ".\n";
            oss << "DOCC and SOCC must specify the occupation of all electrons or none.";
            throw PSIEXCEPTION(oss.str());
        }
    }

    // Check that we have enough basis functions
    for (int h = 0; h < nirrep_; ++h) {
        if (std::max(nalphapi_[h], nbetapi_[h]) > nsopi_[h]) {
            throw PSIEXCEPTION("Not enough basis functions to satisfy requested occupancies");
        }
    }

    // Set additional information
    nuclearrep_ = molecule_->nuclear_repulsion_energy(dipole_field_strength_);
    charge_ = molecule_->molecular_charge();
    multiplicity_ = molecule_->multiplicity();
    nelectron_ = nbeta_ + nalpha_;

    // Copy data for storage
    original_nalphapi_ = nalphapi_;
    original_nbetapi_ = nbetapi_;
    original_nalpha_ = nalpha_;
    original_nbeta_ = nbeta_;

    // Check if it is a broken symmetry solution
    // broken_symmetry_ = false;
    // int socc = 0;
    // for (int h = 0; h < nirrep_; h++) {
    //     socc += soccpi_[h];
    // }

    // if (multiplicity_ == 1 && socc == 2) {
    //     // Set up occupation for the broken symmetry solution
    //     outfile->Printf( "  Broken symmetry solution detected... \n"); //TEST
    //     broken_symmetry_ = true;
    //     int socc_count = 0;
    //     nalphapi_ = doccpi_;
    //     nbetapi_  = doccpi_;

    //     for (int h = 0; h < nirrep_; h++) {
    //         for (int i = 0; i<soccpi_[h]; i++) {
    //             socc_count++;
    //             if (socc_count == 1) {
    //                 nalphapi_[h]++;
    //             }
    //             if (socc_count == 2) {
    //                 nbetapi_[h]++;
    //             }
    //         }
    //     }
    //     if (print_ > 2) {
    //         nalphapi_.print();
    //         nbetapi_.print();
    //     }
    // }

    // How much stuff shall we echo to the user?
    if (options_["PRINT"].has_changed()) print_ = options_.get_int("PRINT");

    initialized_diis_manager_ = false;

    MOM_performed_ = false;  // duplicated py-side (needed before iterate)

    if (print_) {
        print_header();
    }

    // -D is zero by default
    set_scalar_variable("-D Energy", 0.0);  // no-autodoc
    energies_["-D"] = 0.0;

    // CPHF info
    cphf_nfock_builds_ = 0;
    cphf_converged_ = false;
}

void HF::subclass_init() {
    // DFT stuff
    setup_potential();

    if (V_potential() != nullptr) {
        // Do the GRAC
        if (options_.get_double("DFT_GRAC_SHIFT") != 0.0) {
            V_potential()->set_grac_shift(options_.get_double("DFT_GRAC_SHIFT"));
        }

        // Print the KS-specific stuff
        if (print_) {
            V_potential()->print_header();
        }
    }
}

void HF::damping_update(double damping_percentage) {
    throw PSIEXCEPTION(
        "Sorry, damping has not been implemented for this "
        "type of SCF wavefunction yet.");
}

int HF::soscf_update(double soscf_conv, int soscf_min_iter, int soscf_max_iter, int soscf_print) {
    throw PSIEXCEPTION(
        "Sorry, second-order convergence has not been implemented for this "
        "type of SCF wavefunction yet.");
}

void HF::form_V() { throw PSIEXCEPTION("Sorry, DFT functionals are not supported for this type of SCF wavefunction."); }
void HF::form_C(double shift) { throw PSIEXCEPTION("Sorry, the base HF wavefunction cannot construct orbitals."); }
void HF::form_D() { throw PSIEXCEPTION("Sorry, the base HF wavefunction cannot construct densities."); }

std::vector<SharedMatrix> HF::onel_Hx(std::vector<SharedMatrix> x) {
    throw PSIEXCEPTION("Sorry, the base HF wavefunction cannot construct Hx products.");
}
std::vector<SharedMatrix> HF::twoel_Hx(std::vector<SharedMatrix> x, bool combine, std::string return_basis) {
    throw PSIEXCEPTION("Sorry, the base HF wavefunction cannot construct Hx products.");
}
std::vector<SharedMatrix> HF::cphf_Hx(std::vector<SharedMatrix> x) {
    throw PSIEXCEPTION("Sorry, the base HF wavefunction cannot construct cphf_Hx products.");
}
std::vector<SharedMatrix> HF::cphf_solve(std::vector<SharedMatrix> x_vec, double conv_tol, int max_iter,
                                         int print_lvl) {
    throw PSIEXCEPTION("Sorry, the base HF wavefunction cannot solve CPHF equations.");
}
void HF::save_density_and_energy() {
    throw PSIEXCEPTION("Sorry, the base HF wavefunction does not understand a density equation.");
}
void HF::form_G() { throw PSIEXCEPTION("Sorry, the base HF wavefunction does not understand."); }
void HF::form_F() { throw PSIEXCEPTION("Sorry, the base HF wavefunction does not understand Roothan."); }
double HF::compute_E() { throw PSIEXCEPTION("Sorry, the base HF wavefunction does not understand Hall."); }
void HF::rotate_orbitals(SharedMatrix C, const SharedMatrix x) {
    // => Rotate orbitals <= //
    auto U = std::make_shared<Matrix>("Ck", nirrep_, nmopi_, nmopi_);
    std::string reference = options_.get_str("REFERENCE");

    // We guess occ x vir block size by the size of x to make this method easy to use
    Dimension tsize = x->colspi() + x->rowspi();
    if ((reference != "ROHF") && (tsize != nmopi_)) {
        throw PSIEXCEPTION("HF::rotate_orbitals: x dimensions do not match nmo_ dimension.");
    }
    tsize = x->colspi() + x->rowspi() - soccpi();
    if ((reference == "ROHF") && (tsize != nmopi_)) {
        throw PSIEXCEPTION("HF::rotate_orbitals: x dimensions do not match nmo_ dimension.");
    }

    // Form full antisymmetric matrix
    for (size_t h = 0; h < nirrep_; h++) {
        // Whatever the dimension are, we set top right/bot left
        size_t doccpih = (size_t)x->rowspi()[h];
        size_t virpih = (size_t)x->colspi()[h];
        if (!doccpih || !virpih) continue;
        double** up = U->pointer(h);
        double* xp = x->pointer(h)[0];

        // Matrix::schmidt orthogonalizes rows not columns so we need to transpose
        for (size_t i = 0, target = 0; i < doccpih; i++) {
            for (size_t a = (nmopi_[h] - virpih); a < nmopi_[h]; a++) {
                up[a][i] = xp[target];
                up[i][a] = -1.0 * xp[target++];
            }
        }
    }
    U->expm(4, true);

    // Need to build a new one here incase nmo != nso
    auto tmp = linalg::doublet(C, U, false, false);
    C->copy(tmp);
}
void HF::initialize_gtfock_jk() {
    // Build the JK from options, symmetric type
#ifdef HAVE_JK_FACTORY
    // construction of `jk_` depends on communication through legacy_molecule, now removed

    // DGAS is adding to the ghetto, this Python -> C++ -> C -> C++ -> back to C is FUBAR
    if (options_.get_bool("SOSCF"))
        jk_ = std::make_shared<GTFockJK>(basisset_, 2, false);
    else
        jk_ = std::make_shared<GTFockJK>(basisset_, 2, true);
#else
    throw PSIEXCEPTION("GTFock was not compiled in this version.\n");
#endif
}

void HF::finalize() {
    // Clean memory off, handle diis closeout, etc

    // This will be the only one
    if (!options_.get_bool("SAVE_JK")) {
        jk_.reset();
    }

    // Clean up after DIIS
    if (initialized_diis_manager_) diis_manager_.attr("delete_diis_file")();
    diis_manager_ = py::none();
    initialized_diis_manager_ = false;

    // Figure out how many frozen virtual and frozen core per irrep
    compute_fcpi();
    compute_fvpi();
    energy_ = energies_["Total Energy"];

    // Sphalf_.reset();
    X_.reset();
    T_.reset();
}

void HF::set_jk(std::shared_ptr<JK> jk) {
    // Cheap basis check
    int jk_nbf = jk->basisset()->nbf();
    int hf_nbf = basisset_->nbf();
    if (hf_nbf != jk_nbf) {
        throw PSIEXCEPTION("Tried setting a JK object whos number of basis functions does not match HF's!");
    }

    jk_ = jk;
}

void HF::semicanonicalize() { throw PSIEXCEPTION("This type of wavefunction cannot be semicanonicalized!"); }

void HF::find_occupation() {
    // Don't mess with the occ, MOM's got it!
    if (MOM_performed_) {
        MOM();
    } else {
        auto old_nalphapi = nalphapi_;
        auto old_nbetapi = nbetapi_;
        if (!input_docc_ && !input_socc_) {
            assert(nirrep_ == epsilon_a_->nirrep());
            assert(nirrep_ == epsilon_b_->nirrep());

            // The occupations are determined by the Aufbau
            // principle. We first collect all the orbital energies and
            // sort them in increasing order
            std::vector<std::pair<double, int> > pairs_a;
            for (int h = 0; h < nirrep_; ++h) {
                for (int i = 0; i < epsilon_a_->dimpi()[h]; ++i) {
                    pairs_a.push_back(std::make_pair(epsilon_a_->get(h, i), h));
                }
            }
            sort(pairs_a.begin(), pairs_a.end());
            // Same for beta electrons
            std::vector<std::pair<double, int> > pairs_b;
            for (int h = 0; h < nirrep_; ++h) {
                for (int i = 0; i < epsilon_b_->dimpi()[h]; ++i) {
                    pairs_b.push_back(std::make_pair(epsilon_b_->get(h, i), h));
                }
            }
            sort(pairs_b.begin(), pairs_b.end());

            // Sanity check: we must have at least one orbital per electron
            if ((size_t)std::max(nalpha_, nbeta_) > pairs_a.size())
                throw PSIEXCEPTION("Not enough basis functions to satisfy requested occupancies");

            // Reset occupations
            for (int h = 0; h < nirrep_; ++h) {
                nalphapi_[h] = 0;
                nbetapi_[h] = 0;
            }
            // Occupy the lowest nalpha orbitals
            for (int i = 0; i < nalpha_; ++i) nalphapi_[pairs_a[i].second]++;
            // Occupy the lowest nbeta electrons
            for (int i = 0; i < nbeta_; ++i) nbetapi_[pairs_b[i].second]++;
        }

        if (!input_docc_ && !input_socc_) {
            int alphacount = nalphapi_.sum();
            int betacount = nbetapi_.sum();
            if (alphacount != nalpha_) {
                std::ostringstream oss;
                oss << "Count " << alphacount << " alpha electrons, expected " << nalpha_ << ".\n";
                oss << "This is a bug. Please file a report.";
                throw PSIEXCEPTION(oss.str());
            }
            if (betacount != nbeta_) {
                std::ostringstream oss;
                oss << "Count " << betacount << " beta electrons, expected " << nbeta_ << ".\n";
                oss << "This is a bug. Please file a report.";
                throw PSIEXCEPTION(oss.str());
            }
        }

        bool occ_changed = (nalphapi_ != old_nalphapi) || (nbetapi_ != old_nbetapi);

        // If print > 2 (diagnostics), print always
        if ((print_ > 2 || (print_ && occ_changed)) && iteration_ > 0) {
            outfile->Printf("    Occupation by irrep:\n");
            print_occupation();
        }
        // Start MOM if needed (called here because we need the nocc
        // to be decided by Aufbau ordering prior to MOM_start)
        MOM_start();
    }
    // Do fractional orbital normalization here.
    frac();
}

void HF::print_header() {
    int nthread = 1;
#ifdef _OPENMP
    nthread = Process::environment.get_n_threads();
#endif

    outfile->Printf("\n");
    outfile->Printf("         ---------------------------------------------------------\n");
    outfile->Printf("                                   SCF\n");
    outfile->Printf("               by Justin Turney, Rob Parrish, Andy Simmonett\n");
    outfile->Printf("                          and Daniel G. A. Smith\n");
    outfile->Printf("                             %4s Reference\n", options_.get_str("REFERENCE").c_str());
    outfile->Printf("                      %3d Threads, %6ld MiB Core\n", nthread, memory_ / 1048576L);
    outfile->Printf("         ---------------------------------------------------------\n");
    outfile->Printf("\n");
    outfile->Printf("  ==> Geometry <==\n\n");

    molecule_->print();

    outfile->Printf("  Running in %s symmetry.\n\n", molecule_->point_group()->symbol().c_str());

    molecule_->print_rotational_constants();

    outfile->Printf("  Nuclear repulsion = %20.15f\n\n", nuclearrep_);
    outfile->Printf("  Charge       = %d\n", charge_);
    outfile->Printf("  Multiplicity = %d\n", multiplicity_);
    outfile->Printf("  Electrons    = %d\n", nelectron_);
    outfile->Printf("  Nalpha       = %d\n", nalpha_);
    outfile->Printf("  Nbeta        = %d\n\n", nbeta_);

    outfile->Printf("  ==> Algorithm <==\n\n");
    outfile->Printf("  SCF Algorithm Type is %s.\n", options_.get_str("SCF_TYPE").c_str());
    outfile->Printf("  DIIS %s.\n", options_.get_bool("DIIS") ? "enabled" : "disabled");
    if ((options_.get_int("MOM_START") != 0) && (options_["MOM_OCC"].size() != 0))  // TROUBLE, NOT SET YET?
        outfile->Printf("  Excited-state MOM enabled.\n");
    else
        outfile->Printf("  MOM %s.\n", (options_.get_int("MOM_START") == 0) ? "disabled" : "enabled");
    outfile->Printf("  Fractional occupation %s.\n", (options_.get_int("FRAC_START") == 0) ? "disabled" : "enabled");
    outfile->Printf("  Guess Type is %s.\n", options_.get_str("GUESS").c_str());
    outfile->Printf("  Energy threshold   = %3.2e\n", options_.get_double("E_CONVERGENCE"));
    outfile->Printf("  Density threshold  = %3.2e\n", options_.get_double("D_CONVERGENCE"));
    outfile->Printf("  Integral threshold = %3.2e\n\n", integral_threshold_);

    outfile->Printf("  ==> Primary Basis <==\n\n");

    basisset_->print_by_level("outfile", print_);
}

void HF::form_H() {
    T_ = mintshelper()->so_kinetic()->clone();
    V_ = mintshelper()->so_potential()->clone();

    if (debug_ > 2) T_->print("outfile");
    if (debug_ > 2) V_->print("outfile");

    if (perturb_h_) {
        if (dipole_field_type_ == embpot || dipole_field_type_ == sphere ||
            dipole_field_type_ == dx) {  // embedding potential read from file
            if (nirrep_ > 1)
                throw PSIEXCEPTION("RHF_embed: embedding, dx, and spherical potentials require 'symmetry c1'.");
            int nso = 0;
            for (int h = 0; h < nirrep_; h++) nso += nsopi_[h];
            int nao = basisset_->nao();

            // Set up AO->SO transformation matrix (u)
            MintsHelper helper(basisset_, options_, 0);
            SharedMatrix aotoso = helper.petite_list(true)->aotoso();
            Matrix u(nao, nso);
            int offset = 0;

            for (int h = 0; h < nirrep_; h++) {
                // These loops should be vectorized for a (small) efficiency gain.
                for (int j = 0; j < aotoso->coldim(h); j++) {
                    for (int i = 0; i < nao; i++) {
                        u.set(i, j + offset, aotoso->get(h, i, j));
                    }
                }
                offset += aotoso->coldim(h);
            }

            Vector phi_ao(nao);
            Vector phi_so(nso);
            Matrix V_eff(nso, nso);

            if (dipole_field_type_ == embpot) {
                FILE* input = fopen("EMBPOT", "r");
                int npoints;
                int statusvalue = fscanf(input, "%d", &npoints);
                outfile->Printf("  npoints = %d\n", npoints);
                double x, y, z, w, v;
                double max = 0;
                for (int k = 0; k < npoints; k++) {
                    statusvalue = fscanf(input, "%lf %lf %lf %lf %lf", &x, &y, &z, &w, &v);
                    if (std::fabs(v) > max) max = std::fabs(v);

                    basisset_->compute_phi(phi_ao.pointer(), x, y, z);
                    // Transform phi_ao to SO basis
                    phi_so.gemv(true, 1.0, u, phi_ao, 0.0);
                    for (int i = 0; i < nso; i++)
                        for (int j = 0; j < nso; j++) V_eff.add(i, j, w * v * phi_so[i] * phi_so[j]);
                }  // npoints

                outfile->Printf("  Max. embpot value = %20.10f\n", max);
                fclose(input);

            }  // embpot
            else if (dipole_field_type_ == dx) {
                dx_read(V_eff.pointer(), phi_ao.pointer(), phi_so.pointer(), nao, nso, u.pointer());

            }  // dx file
            else if (dipole_field_type_ == sphere) {
                radius_ = options_.get_double("RADIUS");
                thickness_ = options_.get_double("THICKNESS");
                r_points_ = options_.get_int("R_POINTS");
                theta_points_ = options_.get_int("THETA_POINTS");
                phi_points_ = options_.get_int("PHI_POINTS");
                outfile->Printf("  Hard spherical potential radius         = %3.2f bohr\n", radius_);
                outfile->Printf("  Spherical potential thickness           = %3.2f bohr\n", thickness_);
                outfile->Printf("  Number of radial integration points     = %d\n", r_points_);
                outfile->Printf("  Number of colatitude integration points = %d\n", theta_points_);
                outfile->Printf("  Number of azimuthal integration points  = %d\n", phi_points_);

                double r_step = thickness_ / r_points_;         // bohr
                double theta_step = 2 * pc_pi / theta_points_;  // 1 degree in radians
                double phi_step = 2 * pc_pi / phi_points_;      // 1 degree in radians
                double weight = r_step * theta_step * phi_step;
                for (double r = radius_; r < radius_ + thickness_; r += r_step) {
                    for (double theta = 0.0; theta < pc_pi; theta += theta_step) { /* colatitude */
                        for (double phi = 0.0; phi < 2 * pc_pi; phi += phi_step) { /* azimuthal */

                            double x = r * sin(theta) * cos(phi);
                            double y = r * sin(theta) * sin(phi);
                            double z = r * cos(theta);

                            double jacobian = weight * r * r * sin(theta);

                            basisset_->compute_phi(phi_ao.pointer(), x, y, z);

                            phi_so.gemv(true, 1.0, u, phi_ao, 0.0);

                            for (int i = 0; i < nso; i++)
                                for (int j = 0; j < nso; j++)
                                    V_eff.add(i, j, jacobian * (-1.0e6) * phi_so[i] * phi_so[j]);
                        }
                    }
                }
            }  // sphere

            outfile->Printf("  Perturbing H by %f %f %f V_eff.\n", dipole_field_strength_[0], dipole_field_strength_[1],
                            dipole_field_strength_[2]);
            if (options_.get_int("PRINT") > 3) V_eff.print_out();

            if (dipole_field_type_ == dx) {
                V_->copy(V_eff);
            } else {
                V_->add(V_eff);
            }

        }  // embpot or sphere
    }      // end perturb_h_

    // If an external field exists, add it to the one-electron Hamiltonian
    if (external_pot_) {
        if (options_.get_bool("EXTERNAL_POTENTIAL_SYMMETRY") == false && H_->nirrep() != 1)
            throw PSIEXCEPTION("SCF: External Fields are not consistent with symmetry. Set symmetry c1.");

        auto Vprime = external_pot_->computePotentialMatrix(basisset_);

        if (options_.get_bool("EXTERNAL_POTENTIAL_SYMMETRY")) {
            // Attempt to apply symmetry. No error checking is performed.
            auto Vprimesym = factory_->create_shared_matrix("External Potential");
            Vprimesym->apply_symmetry(Vprime, AO2SO_);
            Vprime = Vprimesym;
        }

        if (print_) {
            external_pot_->set_print(print_);
            external_pot_->print();
        }
        if (print_ > 3) Vprime->print();
        V_->add(Vprime);

        // Extra nuclear repulsion
        double enuc2 = external_pot_->computeNuclearEnergy(molecule_);
        if (print_) {
            outfile->Printf("  Old nuclear repulsion        = %20.15f\n", nuclearrep_);
            outfile->Printf("  Additional nuclear repulsion = %20.15f\n", enuc2);
            outfile->Printf("  Total nuclear repulsion      = %20.15f\n\n", nuclearrep_ + enuc2);
        }
        nuclearrep_ += enuc2;

    }  // end external

    // Save perturbed V_ for future (e.g. correlated) calcs
    V_->save(psio_, PSIF_OEI);

    H_->copy(T_);
    H_->add(V_);

    if (print_ > 3) H_->print("outfile");
}

void HF::form_Shalf() {
    BasisSetOrthogonalization::OrthogonalizationMethod method;
    if (options_.get_str("S_ORTHOGONALIZATION") == "SYMMETRIC")
        method = BasisSetOrthogonalization::Symmetric;
    else if (options_.get_str("S_ORTHOGONALIZATION") == "CANONICAL")
        method = BasisSetOrthogonalization::Canonical;
    else if (options_.get_str("S_ORTHOGONALIZATION") == "PARTIALCHOLESKY")
        method = BasisSetOrthogonalization::PartialCholesky;
    else if (options_.get_str("S_ORTHOGONALIZATION") == "AUTO")
        method = BasisSetOrthogonalization::Automatic;
    else
        throw PSIEXCEPTION("Unrecognized S_ORTHOGONALIZATION method\n");

    bool used_brian = false;

#if USING_BrianQC
    if (brianEnable) {
        double S_cutoff = options_.get_double("S_TOLERANCE");
        if (print_)
            outfile->Printf("  BrianQC enabled, using Canonical Orthogonalization with cutoff of %14.10E.\n", S_cutoff);

        brianInt computeOverlapRoot = BRIAN_FALSE;
        brianInt computeOverlapInverseRoot = BRIAN_TRUE;
        brianInt basisRank;
        SharedMatrix buffer = std::make_shared<Matrix>(nirrep_, nsopi_, nsopi_);
        brianSCFComputeOverlapRoot(&brianCookie, &computeOverlapRoot, &computeOverlapInverseRoot, S_->get_pointer(),
                                   &S_cutoff, &basisRank, nullptr, buffer->get_pointer());
        checkBrian();

        nmo_ = basisRank;
        nmopi_[0] = basisRank;

        X_->init(nirrep_, nsopi_, nmopi_, "X (Canonical Orthogonalization)");
        for (int i = 0; i < nso_; i++) {
            for (int j = 0; j < nmo_; j++) {
                X_->set(i, j, buffer->get(nmo_ - 1 - j, i));
            }
        }

        if (print_) outfile->Printf("  Overall, %d of %d possible MOs eliminated.\n\n", nso_ - nmo_, nso_);

        used_brian = true;
    }
#endif

    if (!used_brian) {
        double lindep_tolerance = options_.get_double("S_TOLERANCE");
        double cholesky_tolerance = options_.get_double("S_CHOLESKY_TOLERANCE");

        BasisSetOrthogonalization orthog(method, S_, lindep_tolerance, cholesky_tolerance, print_);

        // Transform
        X_ = orthog.basis_to_orthog_basis();

        // Update nmo_
        nmopi_ = X_->colspi();
        nmo_ = nmopi_.sum();
    }

    // Double check occupation vectors
    for (int h = 0; h < X_->nirrep(); ++h) {
        if (std::max(nalphapi_[h], nbetapi_[h]) > nmopi_[h]) {
            throw PSIEXCEPTION("Not enough molecular orbitals to satisfy requested occupancies");
        }
    }
    // Refreshes twice in RHF, no big deal
    epsilon_a_->init(nmopi_);
    Ca_->init(nirrep_, nsopi_, nmopi_, "Alpha MO coefficients");
    epsilon_b_->init(nmopi_);
    if (!same_a_b_orbs_) {
        Cb_->init(nirrep_, nsopi_, nmopi_, "Beta MO coefficients");
    }

    // Extra matrix dimension changes for specific derived classes
    prepare_canonical_orthogonalization();

    if (print_ > 3) {
        S_->print("outfile");
        X_->print("outfile");
    }
}

void HF::compute_fcpi() {
    // FROZEN_DOCC takes precedence, FREEZE_CORE directive has second priority
    if (options_["FROZEN_DOCC"].has_changed()) {
        if (options_["FROZEN_DOCC"].size() != epsilon_a_->nirrep()) {
            throw PSIEXCEPTION("The FROZEN_DOCC array has the wrong dimensions");
        }
        for (int h = 0; h < epsilon_a_->nirrep(); h++) {
            frzcpi_[h] = options_["FROZEN_DOCC"][h].to_integer();
        }
    } else {
        int nfzc = 0;
        if (options_.get_int("NUM_FROZEN_DOCC") != 0) {
            nfzc = options_.get_int("NUM_FROZEN_DOCC");
        } else {
            nfzc = basisset_->n_frozen_core();
        }
        // Print out orbital energies.
        std::vector<std::pair<double, int> > pairs;
        for (int h = 0; h < epsilon_a_->nirrep(); ++h) {
            for (int i = 0; i < epsilon_a_->dimpi()[h]; ++i) pairs.push_back(std::make_pair(epsilon_a_->get(h, i), h));
            frzcpi_[h] = 0;
        }
        sort(pairs.begin(), pairs.end());

        for (int i = 0; i < nfzc; ++i) frzcpi_[pairs[i].second]++;
    }
    // total frozen core
    nfrzc_ = 0;
    for (int h = 0; h < epsilon_a_->nirrep(); h++) nfrzc_ += frzcpi_[h];
}

void HF::compute_fvpi() {
    // FROZEN_UOCC takes precedence, FREEZE_UOCC directive has second priority
    if (options_["FROZEN_UOCC"].has_changed()) {
        if (options_["FROZEN_UOCC"].size() != epsilon_a_->nirrep()) {
            throw PSIEXCEPTION("The FROZEN_UOCC array has the wrong dimensions");
        }
        for (int h = 0; h < epsilon_a_->nirrep(); h++) {
            frzvpi_[h] = options_["FROZEN_UOCC"][h].to_integer();
        }
    } else {
        int nfzv = options_.get_int("NUM_FROZEN_UOCC");
        // Print out orbital energies.
        std::vector<std::pair<double, int> > pairs;
        for (int h = 0; h < epsilon_a_->nirrep(); ++h) {
            for (int i = 0; i < epsilon_a_->dimpi()[h]; ++i) pairs.push_back(std::make_pair(epsilon_a_->get(h, i), h));
            frzvpi_[h] = 0;
        }
        sort(pairs.begin(), pairs.end(), std::greater<std::pair<double, int> >());

        for (int i = 0; i < nfzv; ++i) frzvpi_[pairs[i].second]++;
    }
}

void HF::print_orbital_pairs(const char* header, std::vector<std::pair<double, std::pair<std::string, int> > > orbs) {
    outfile->Printf("    %-70s\n\n    ", header);
    int count = 0;
    for (int i = 0; i < orbs.size(); i++) {
        outfile->Printf("%4d%-4s%11.6f  ", orbs[i].second.second, orbs[i].second.first.c_str(), orbs[i].first);
        if (count++ % 3 == 2 && count != orbs.size()) outfile->Printf("\n    ");
    }
    outfile->Printf("\n\n");
}

void HF::print_orbitals() {
    std::vector<std::string> labels = molecule_->irrep_labels();

    outfile->Printf("    Orbital Energies [Eh]\n    ---------------------\n\n");

    std::string reference = options_.get_str("REFERENCE");
    if ((reference == "RHF") || (reference == "RKS")) {
        std::vector<std::pair<double, std::pair<std::string, int> > > occ;
        std::vector<std::pair<double, std::pair<std::string, int> > > vir;

        for (int h = 0; h < nirrep_; h++) {
            std::vector<std::pair<double, int> > orb_e;
            for (int a = 0; a < nmopi_[h]; a++) orb_e.push_back(std::make_pair(epsilon_a_->get(h, a), a));
            std::sort(orb_e.begin(), orb_e.end());

            std::vector<int> orb_order(nmopi_[h]);
            for (int a = 0; a < nmopi_[h]; a++) orb_order[orb_e[a].second] = a;

            for (int a = 0; a < nalphapi_[h]; a++)
                occ.push_back(std::make_pair(epsilon_a_->get(h, a), std::make_pair(labels[h], orb_order[a] + 1)));
            for (int a = nalphapi_[h]; a < nmopi_[h]; a++)
                vir.push_back(std::make_pair(epsilon_a_->get(h, a), std::make_pair(labels[h], orb_order[a] + 1)));
        }
        std::sort(occ.begin(), occ.end());
        std::sort(vir.begin(), vir.end());

        print_orbital_pairs("Doubly Occupied:", occ);
        print_orbital_pairs("Virtual:", vir);

    } else if ((reference == "UHF") || (reference == "UKS") || (reference == "CUHF")) {
        std::vector<std::pair<double, std::pair<std::string, int> > > occA;
        std::vector<std::pair<double, std::pair<std::string, int> > > virA;
        std::vector<std::pair<double, std::pair<std::string, int> > > occB;
        std::vector<std::pair<double, std::pair<std::string, int> > > virB;

        for (int h = 0; h < nirrep_; h++) {
            std::vector<std::pair<double, int> > orb_eA;
            for (int a = 0; a < nmopi_[h]; a++) orb_eA.push_back(std::make_pair(epsilon_a_->get(h, a), a));
            std::sort(orb_eA.begin(), orb_eA.end());

            std::vector<int> orb_orderA(nmopi_[h]);
            for (int a = 0; a < nmopi_[h]; a++) orb_orderA[orb_eA[a].second] = a;

            for (int a = 0; a < nalphapi_[h]; a++)
                occA.push_back(std::make_pair(epsilon_a_->get(h, a), std::make_pair(labels[h], orb_orderA[a] + 1)));
            for (int a = nalphapi_[h]; a < nmopi_[h]; a++)
                virA.push_back(std::make_pair(epsilon_a_->get(h, a), std::make_pair(labels[h], orb_orderA[a] + 1)));

            std::vector<std::pair<double, int> > orb_eB;
            for (int a = 0; a < nmopi_[h]; a++) orb_eB.push_back(std::make_pair(epsilon_b_->get(h, a), a));
            std::sort(orb_eB.begin(), orb_eB.end());

            std::vector<int> orb_orderB(nmopi_[h]);
            for (int a = 0; a < nmopi_[h]; a++) orb_orderB[orb_eB[a].second] = a;

            for (int a = 0; a < nbetapi_[h]; a++)
                occB.push_back(std::make_pair(epsilon_b_->get(h, a), std::make_pair(labels[h], orb_orderB[a] + 1)));
            for (int a = nbetapi_[h]; a < nmopi_[h]; a++)
                virB.push_back(std::make_pair(epsilon_b_->get(h, a), std::make_pair(labels[h], orb_orderB[a] + 1)));
        }
        std::sort(occA.begin(), occA.end());
        std::sort(virA.begin(), virA.end());
        std::sort(occB.begin(), occB.end());
        std::sort(virB.begin(), virB.end());

        print_orbital_pairs("Alpha Occupied:", occA);
        print_orbital_pairs("Alpha Virtual:", virA);
        print_orbital_pairs("Beta Occupied:", occB);
        print_orbital_pairs("Beta Virtual:", virB);

    } else if (reference == "ROHF") {
        std::vector<std::pair<double, std::pair<std::string, int> > > docc;
        std::vector<std::pair<double, std::pair<std::string, int> > > socc;
        std::vector<std::pair<double, std::pair<std::string, int> > > vir;

        for (int h = 0; h < nirrep_; h++) {
            std::vector<std::pair<double, int> > orb_e;
            for (int a = 0; a < nmopi_[h]; a++) orb_e.push_back(std::make_pair(epsilon_a_->get(h, a), a));
            std::sort(orb_e.begin(), orb_e.end());

            std::vector<int> orb_order(nmopi_[h]);
            for (int a = 0; a < nmopi_[h]; a++) orb_order[orb_e[a].second] = a;

            for (int a = 0; a < nbetapi_[h]; a++)
                docc.push_back(std::make_pair(epsilon_a_->get(h, a), std::make_pair(labels[h], orb_order[a] + 1)));
            for (int a = nbetapi_[h]; a < nalphapi_[h]; a++)
                socc.push_back(std::make_pair(epsilon_a_->get(h, a), std::make_pair(labels[h], orb_order[a] + 1)));
            for (int a = nalphapi_[h]; a < nmopi_[h]; a++)
                vir.push_back(std::make_pair(epsilon_a_->get(h, a), std::make_pair(labels[h], orb_order[a] + 1)));
        }
        std::sort(docc.begin(), docc.end());
        std::sort(socc.begin(), socc.end());
        std::sort(vir.begin(), vir.end());

        print_orbital_pairs("Doubly Occupied:", docc);
        print_orbital_pairs("Singly Occupied:", socc);
        print_orbital_pairs("Virtual:", vir);

    } else {
        throw PSIEXCEPTION("Unknown reference in HF::print_orbitals");
    }

    outfile->Printf("    Final Occupation by Irrep:\n");
    print_occupation();
}

void HF::compute_sapgau_guess() {
  // Build auxiliary basis set object
  auto sap_basis = get_basisset("SAPGAU");
  // Do the SAP magic to the basis
  sap_basis->convert_sap_contraction();

  auto zero_basis = BasisSet::zero_ao_basis_set();
  auto nsap = sap_basis->nbf();
  auto nbf = basisset_->nbf();

  // Build (P|pq) raw 3-index ERIs in AO basis, dimension (Nsap, 1, nbf, nbf).
  auto Ppq = mintshelper()->ao_eri(sap_basis, zero_basis, basisset_, basisset_);

  // Build repulsive potential matrix in AO basis.
  auto Vsap = std::make_shared<Matrix>("VSAP", basisset_->nbf(), basisset_->nbf());
  auto Varr = Vsap->pointer();
  auto Parr = Ppq->pointer();
  for (auto P = 0; P < nsap; P++) {
    for(auto u = 0; u < nbf; u++) {
      for(auto v = 0; v < nbf; v++) {
        // TBD: this is using Natoms times too much memory - the
        // integrals should be computed in a loop, parallellizing over
        // the AO shell pairs
        Varr[u][v] += Parr[P][u*nbf+v];
      }
    }
  }

  // Convert repulsive potential into the SO basis
  auto Fsap = std::make_shared<Matrix>("FSAP", AO2SO_->colspi(), AO2SO_->colspi());
  Fsap->apply_symmetry(Vsap, AO2SO_);
  // and add in the core Hamiltonian
  Fsap->add(H_);

  // Set the alpha and beta Fock matrices
  Fa_->copy(Fsap);
  Fb_->copy(Fsap);
}

void HF::guess() {
    // don't save guess energy as "the" energy because we need to avoid
    // a false positive test for convergence on the first iteration (that
    // was happening before in tests/scf-guess-read before I removed
    // the statements putting this into E_).  -CDS 3/25/13
    double guess_E;

    // What does the user want?
    // Options will be:
    // "CORE"-CORE Hamiltonain
    // "GWH"-Generalized Wolfsberg-Helmholtz
    // "SAD"-Superposition of Atomic Densities
    std::string guess_type = options_.get_str("GUESS");

    // Take care of options that should be overridden
    if (guess_type == "AUTO") {
        outfile->Printf("\nWarning! Guess was AUTO, switching to CORE!\n\n");
        outfile->Printf("           This option should have been configured at the driver level.\n\n");
        guess_type = "CORE";
    }

    if ((guess_type == "READ") && !guess_Ca_) {
        outfile->Printf("\nWarning! Guess was READ without Ca set, switching to CORE!\n");
        outfile->Printf("           This option should have been configured at the driver level.\n\n");
        guess_type = "CORE";
    }

    if (guess_Ca_) {
        if (print_) outfile->Printf("  SCF Guess: Orbitals guess was supplied from a previous computation.\n\n");

        std::string reference = options_.get_str("REFERENCE");
        bool single_orb = (reference == "RHF");

        if (single_orb) {
            guess_Cb_ = guess_Ca_;
        } else {
            if (!guess_Cb_) {
                throw PSIEXCEPTION("Guess Ca was set, but did not find a matching Cb!\n");
            }
        }

        if ((guess_Ca_->nirrep() != nirrep_) || (guess_Cb_->nirrep() != nirrep_)) {
            throw PSIEXCEPTION(
                "Number of guess of the input orbitals do not match number of irreps of the wavefunction.");
        }
        if ((guess_Ca_->rowspi() != nsopi_) || (guess_Cb_->rowspi() != nsopi_)) {
            throw PSIEXCEPTION("Nso of the guess orbitals do not match Nso of the wavefunction.");
        }

        for (int h = 0; h < nirrep_; h++) {
            for (int i = 0; i < guess_Ca_->colspi()[h]; i++) {
                C_DCOPY(nsopi_[h], &guess_Ca_->pointer(h)[0][i], guess_Ca_->colspi()[h], &Ca_->pointer(h)[0][i],
                        nmopi_[h]);
            }
        }

        if (single_orb) {
            Cb_ = Ca_;
        } else {
            for (int h = 0; h < nirrep_; h++) {
                for (int i = 0; i < guess_Cb_->colspi()[h]; i++) {
                    C_DCOPY(nsopi_[h], &guess_Cb_->pointer(h)[0][i], guess_Cb_->colspi()[h], &Cb_->pointer(h)[0][i],
                            nmopi_[h]);
                }
            }
        }

        // Figure out occupations from given input
        if (!(input_socc_ || input_docc_)) {
            nalphapi_ = guess_Ca_->colspi();
            nbetapi_ = guess_Cb_->colspi();
            nalpha_ = nalphapi_.sum();
            nbeta_ = nbetapi_.sum();
        }

        format_guess();
        form_D();

        // This is a guess iteration: orbital occupations may be reset in SCF
        iteration_ = -1;
        guess_E = compute_initial_E();

    } else if (guess_type == "SAD") {
        if (print_)
            outfile->Printf(
                "  SCF Guess: Superposition of Atomic Densities via on-the-fly atomic UHF (no occupation "
                "information).\n\n");

        // Superposition of Atomic Density. Modified by Susi Lehtola
        // 2018-12-15 to work also for ROHF, as well as to allow using
        // SAD with predefined orbital occupations. The algorithm is
        // the same as in van Lenthe et al, "Starting SCF Calculations
        // by Superposition of Atomic Densities", J Comput Chem 27,
        // 926 (2006).

        // Build non-idempotent, spin-restricted SAD density matrix
        compute_SAD_guess(false);

        // This is a guess iteration: orbital occupations must be
        // reset in SCF.
        iteration_ = -1;
        // SAD doesn't yield orbitals so also the SCF logic is
        // slightly different for the first iteration.
        sad_ = true;
        guess_E = compute_initial_E();

    } else if (guess_type == "SADNO") {
        if (print_)
            outfile->Printf(
                "  SCF Guess: Superposition of Atomic Densities' Natural Orbitals via on-the-fly atomic UHF "
                "(doi:10.1021/acs.jctc.8b01089).\n\n");

        // Like the above, but builds natural orbitals from the SAD
        // density matrix.

        // Build non-idempotent, spin-restricted SAD density matrix
        compute_SAD_guess(true);
        // Find occupations
        find_occupation();

        // Now we have orbitals and occupations, build a density matrix
        form_D();
        guess_E = compute_initial_E();

    } else if (guess_type == "HUCKEL") {
        if (print_)
            outfile->Printf("  SCF Guess: Huckel guess via on-the-fly atomic UHF (doi:10.1021/acs.jctc.8b01089).\n\n");

        // Huckel guess, written by Susi Lehtola 2019-01-27.  See "An
        // assessment of initial guesses for self-consistent field
        // calculations. Superposition of Atomic Potentials: simple
        // yet efficient", JCTC 2019, doi: 10.1021/acs.jctc.8b01089.

        if (!options_.get_bool("SAD_SPIN_AVERAGE")) {
            throw PSIEXCEPTION("  Huckel guess requires SAD_SPIN_AVERAGE = True!");
        }
        if (!options_.get_bool("SAD_FRAC_OCC")) {
            throw PSIEXCEPTION("  Huckel guess requires SAD_FRAC_OCC = True!");
        }
        compute_huckel_guess(false);

        form_initial_C();
        form_D();
        guess_E = compute_initial_E();

    } else if (guess_type == "MODHUCKEL") {
      if (print_)
            outfile->Printf("  SCF Guess: Huckel guess via on-the-fly atomic UHF (doi:10.1021/acs.jctc.8b01089) with the updated GWH rule from doi:10.1021/ja00480a005.\n\n");

        // Huckel guess, written by Susi Lehtola 2019-01-27.  See "An
        // assessment of initial guesses for self-consistent field
        // calculations. Superposition of Atomic Potentials: simple
        // yet efficient", JCTC 2019, doi: 10.1021/acs.jctc.8b01089.

        if (!options_.get_bool("SAD_SPIN_AVERAGE")) {
            throw PSIEXCEPTION("  Huckel guess requires SAD_SPIN_AVERAGE = True!");
        }
        if (!options_.get_bool("SAD_FRAC_OCC")) {
            throw PSIEXCEPTION("  Huckel guess requires SAD_FRAC_OCC = True!");
        }
        compute_huckel_guess(true);

        form_initial_C();
        form_D();
        guess_E = compute_initial_E();

    } else if (guess_type == "GWH") {
        // Generalized Wolfsberg Helmholtz (Sounds cool, easy to code)
        if (print_) outfile->Printf("  SCF Guess: Generalized Wolfsberg-Helmholtz applied to core Hamiltonian.\n\n");

        Fa_->zero();  // Try Fa_{mn} = S_{mn} (H_{mm} + H_{nn})/2
        int h, i, j;
        const int* opi = S_->rowspi();
        int nirreps = S_->nirrep();
        for (h = 0; h < nirreps; ++h) {
            for (i = 0; i < opi[h]; ++i) {
                Fa_->set(h, i, i, H_->get(h, i, i));
                for (j = 0; j < i; ++j) {
                    Fa_->set(h, i, j, 0.875 * S_->get(h, i, j) * (H_->get(h, i, i) + H_->get(h, j, j)));
                    Fa_->set(h, j, i, Fa_->get(h, i, j));
                }
            }
        }
        Fb_->copy(Fa_);
        form_initial_C();
        form_D();
        guess_E = compute_initial_E();

    } else if (guess_type == "CORE") {
        if (print_) outfile->Printf("  SCF Guess: Core (One-Electron) Hamiltonian.\n\n");

        Fa_->copy(H_);  // Try the core Hamiltonian as the Fock Matrix
        Fb_->copy(H_);

        form_initial_C();
        form_D();
        guess_E = compute_initial_E();

    } else if (guess_type == "SAP") {
        // SAP guess
        if (print_)
            outfile->Printf("  SCF Guess: Superposition of Atomic Potentials (doi:10.1021/acs.jctc.8b01089).\n\n");

        auto builder = VBase::build_V(basisset_, functional_, options_, "SAP");
        builder->initialize();

        // Print info on the integration grid
        if (print_) {
            builder->print_header();
        }

        // Build the SAP potential
        std::vector<SharedMatrix> Vsap;
        Vsap.push_back(SharedMatrix(factory_->create_matrix("Vsap")));
        builder->compute_V(Vsap);
        Fa_->copy(T_);
        Fa_->add(Vsap[0]);
        Fb_->copy(Fa_);
        form_initial_C();
        form_D();
        guess_E = compute_initial_E();

    } else if (guess_type == "SAPGAU") {
      if (print_)
        outfile->Printf("  SCF Guess: Superposition of Atomic Potentials (doi:10.1021/acs.jctc.8b01089).\n  Using error function fits of the atomic potentials (doi:10.1063/5.0004046).\n\n");

        // Build the SAP potential
        compute_sapgau_guess();
        form_initial_C();

        // Find occupations
        find_occupation();

        // Now we have orbitals and occupations, build a density matrix
        form_D();
        guess_E = compute_initial_E();

    } else {
        throw PSIEXCEPTION("  SCF Guess: No guess was found!");
    }

    if (print_ > 3) {
        Ca_->print();
        Cb_->print();
        Da_->print();
        Db_->print();
        Fa_->print();
        Fb_->print();
    }
    energies_["Total Energy"] = 0.0;  // don't use this guess in our convergence checks
}

void HF::format_guess() {
    // Nothing to do, only for special cases
}

void HF::check_phases() {
    for (int h = 0; h < nirrep_; ++h) {
        for (int p = 0; p < Ca_->colspi(h); ++p) {
            for (int mu = 0; mu < Ca_->rowspi(h); ++mu) {
                if (std::fabs(Ca_->get(h, mu, p)) > 1.0E-3) {
                    if (Ca_->get(h, mu, p) < 1.0E-3) {
                        Ca_->scale_column(h, p, -1.0);
                    }
                    break;
                }
            }
        }
    }

    if (Ca_ != Cb_) {
        for (int h = 0; h < nirrep_; ++h) {
            for (int p = 0; p < Cb_->colspi(h); ++p) {
                for (int mu = 0; mu < Cb_->rowspi(h); ++mu) {
                    if (std::fabs(Cb_->get(h, mu, p)) > 1.0E-3) {
                        if (Cb_->get(h, mu, p) < 1.0E-3) {
                            Cb_->scale_column(h, p, -1.0);
                        }
                        break;
                    }
                }
            }
        }
    }
}

void HF::print_occupation() {
    auto labels = molecule_->irrep_labels();
    auto reference = options_.get_str("REFERENCE");
    outfile->Printf("          ");
    for (int h = 0; h < nirrep_; ++h) outfile->Printf(" %4s ", labels[h].c_str());
    outfile->Printf("\n");
    auto docc = doccpi();
    auto socc = soccpi();
    outfile->Printf("    DOCC [ ");
    for (int h = 0; h < nirrep_ - 1; ++h) outfile->Printf(" %4d,", docc[h]);
    outfile->Printf(" %4d ]\n", docc[nirrep_ - 1]);
    if (reference != "RHF" && reference != "RKS") {
        outfile->Printf("    SOCC [ ");
        for (int h = 0; h < nirrep_ - 1; ++h) outfile->Printf(" %4d,", socc[h]);
        outfile->Printf(" %4d ]\n", socc[nirrep_ - 1]);
    }
    // Also print nalpha and nbeta per irrep, which are more physically meaningful
    outfile->Printf("    NA   [ ");
    for (int h = 0; h < nirrep_ - 1; ++h) outfile->Printf(" %4d,", nalphapi_[h]);
    outfile->Printf(" %4d ]\n", nalphapi_[nirrep_ - 1]);
    outfile->Printf("    NB   [ ");
    for (int h = 0; h < nirrep_ - 1; ++h) outfile->Printf(" %4d,", nbetapi_[h]);
    outfile->Printf(" %4d ]\n", nbetapi_[nirrep_ - 1]);

    outfile->Printf("\n");
}

//  Returns a vector of the occupation of the a orbitals
std::shared_ptr<Vector> HF::occupation_a() const {
    auto occA = std::make_shared<Vector>(nmopi_);
    for (int h = 0; h < nirrep_; ++h)
        for (int n = 0; n < nalphapi()[h]; n++) occA->set(h, n, 1.0);

    return occA;
}

//  Returns a vector of the occupation of the b orbitals
std::shared_ptr<Vector> HF::occupation_b() const {
    auto occB = std::make_shared<Vector>(nmopi_);
    for (int h = 0; h < nirrep_; ++h)
        for (int n = 0; n < nbetapi()[h]; n++) occB->set(h, n, 1.0);

    return occB;
}

void HF::diagonalize_F(const SharedMatrix& Fm, SharedMatrix& Cm, std::shared_ptr<Vector>& epsm) {
#ifdef USING_BrianQC
    if (brianEnable) {
        brianInt basisSize = basisset_->nbf();
        brianInt basisRank = X_->coldim(0);

        // BrianQC needs the matrices in a column-major memory layout,
        // so we construct a temporary transposed version of the X matrix,
        // and allocate a transposed C matrix for BrianQC to fill
        std::shared_ptr<Matrix> orthonormalizationMatrix = X_->transpose();
        std::shared_ptr<Matrix> C = std::make_shared<Matrix>(basisRank, basisSize);

        brianSCFDiagonalizeFock(&brianCookie, &basisRank, Fm->get_pointer(0), orthonormalizationMatrix->get_pointer(0),
                                C->get_pointer(0), epsm->pointer(0));
        checkBrian();

        Cm->copy(C->transpose());

        return;
    }
#endif

    // Form F' = X'FX for canonical orthogonalization
    auto diag_F_temp = linalg::triplet(X_, Fm, X_, true, false, false);

    // Form C' = eig(F')
    auto diag_C_temp = std::make_shared<Matrix>(nirrep_, nmopi_, nmopi_);
    diag_F_temp->diagonalize(diag_C_temp, epsm);

    // Form C = XC'
    Cm->gemm(false, false, 1.0, X_, diag_C_temp, 0.0);
}

void HF::reset_occupation() {
    // RHF style for now
    nalphapi_ = original_nalphapi_;
    nbetapi_ = original_nbetapi_;

    // These may not match the per irrep. Will remap correctly next find_occupation call
    nalpha_ = original_nalpha_;
    nbeta_ = original_nbeta_;
}

SharedMatrix HF::form_Fia(SharedMatrix Fso, SharedMatrix Cso, int* noccpi) {
    const int* nsopi = Cso->rowspi();
    const int* nmopi = Cso->colspi();
    int* nvirpi = new int[nirrep_];

    for (int h = 0; h < nirrep_; h++) nvirpi[h] = nmopi[h] - noccpi[h];

    auto Fia = std::make_shared<Matrix>("Fia (Some Basis)", nirrep_, noccpi, nvirpi);

    // Hack to get orbital e for this Fock
    auto C2 = std::make_shared<Matrix>("C2", Cso->rowspi(), Cso->colspi());
    auto E2 = std::make_shared<Vector>("E2", Cso->colspi());
    diagonalize_F(Fso, C2, E2);

    for (int h = 0; h < nirrep_; h++) {
        int nmo = nmopi[h];
        int nso = nsopi[h];
        int nvir = nvirpi[h];
        int nocc = noccpi[h];

        if (nmo == 0 || nso == 0 || nvir == 0 || nocc == 0) continue;

        // double** C = Cso->pointer(h);
        double** C = C2->pointer(h);
        double** F = Fso->pointer(h);
        double** Fiap = Fia->pointer(h);

        double** Temp = block_matrix(nocc, nso);

        C_DGEMM('T', 'N', nocc, nso, nso, 1.0, C[0], nmo, F[0], nso, 0.0, Temp[0], nso);
        C_DGEMM('N', 'N', nocc, nvir, nso, 1.0, Temp[0], nso, &C[0][nocc], nmo, 0.0, Fiap[0], nvir);

        free_block(Temp);

        // double* eps = E2->pointer(h);
        // for (int i = 0; i < nocc; i++)
        //    for (int a = 0; a < nvir; a++)
        //        Fiap[i][a] /= eps[a + nocc] - eps[i];
    }

    // Fia->print();

    delete[] nvirpi;

    return Fia;
}
SharedMatrix HF::form_FDSmSDF(SharedMatrix Fso, SharedMatrix Dso) {
    auto FDSmSDF = linalg::triplet(Fso, Dso, S_, false, false, false);
    auto SDF = FDSmSDF->transpose();
    FDSmSDF->subtract(SDF);

    SDF.reset();

    FDSmSDF->transform(X_);

    return FDSmSDF;
}

void HF::print_stability_analysis(std::vector<std::pair<double, int> >& vec) const {
    std::sort(vec.begin(), vec.end());
    std::vector<std::pair<double, int> >::const_iterator iter = vec.begin();
    outfile->Printf("    ");
    auto irrep_labels = molecule_->irrep_labels();
    int count = 0;
    for (; iter != vec.end(); ++iter) {
        ++count;
        outfile->Printf("%4s %-10.6f", irrep_labels[iter->second].c_str(), iter->first);
        if (count == 4) {
            outfile->Printf("\n    ");
            count = 0;
        } else {
            outfile->Printf("    ");
        }
    }
    if (count)
        outfile->Printf("\n\n");
    else
        outfile->Printf("\n");
}
bool HF::stability_analysis() {
    throw PSIEXCEPTION("Stability analysis hasn't been implemented yet for this wfn type.");
    return false;
}
}  // namespace scf
}  // namespace psi
