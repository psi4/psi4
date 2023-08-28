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
// Need libint for maximum angular momentum
#ifdef ENABLE_Libint1t
#include <libint/libint.h>
#endif
#include <libint2/shell.h>
/*!
    \defgroup MINTS libmints: Integral library
    \ingroup MINTS
*/

#include "psi4/libciomr/libciomr.h"
#include "psi4/psifiles.h"
#include "psi4/libpsi4util/libpsi4util.h"

#include "vector3.h"
#include "molecule.h"
#include "basisset.h"
#include "dimension.h"
#include "sobasis.h"
#include "integral.h"
#include "gshell.h"
#include "factory.h"
#include "pointgrp.h"
#include "wavefunction.h"
#include "coordentry.h"
#include "psi4/libpsi4util/process.h"

#include <memory>
#include <regex>
#include <stdexcept>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <map>
#include <list>

using namespace psi;

bool BasisSet::initialized_shared_ = false;

std::vector<Vector3> BasisSet::exp_ao[LIBINT_MAX_AM];

namespace {
bool has_ending(std::string const &fullString, std::string const &ending) {
    if (fullString.length() >= ending.length()) {
        return (0 == fullString.compare(fullString.length() - ending.length(), ending.length(), ending));
    } else {
        return false;
    }
}

std::string to_upper_copy(const std::string &original) {
    std::string upper = original;
    to_upper(upper);
    return upper;
}
}  // namespace

// Constructs a zero AO basis set
BasisSet::BasisSet() {
    initialize_singletons();

    // Add a dummy atom at the origin, to hold this basis function
    molecule_ = std::make_shared<Molecule>();
    molecule_->add_atom(0, 0.0, 0.0, 0.0);
    // Fill with data representing a single S function, at the origin, with 0 exponent
    n_uprimitive_ = 1;
    n_shells_ = 1;
    nprimitive_ = 1;
    nao_ = 1;
    nbf_ = 1;
    n_prim_per_shell_ = std::vector<int>(1, 1);
    uexponents_ = std::vector<double>(1, 0.0);
    ucoefficients_ = std::vector<double>(1, 1.0);
    uerd_coefficients_ = std::vector<double>(1, 1.0);
    uoriginal_coefficients_ = std::vector<double>(1, 1.0);
    shell_first_ao_ = std::vector<int>(1, 0);
    shell_first_basis_function_ = std::vector<int>(1, 0);
    shells_ = std::vector<GaussianShell>(1);
    l2_shells_.push_back(libint2::Shell::unit());
    ao_to_shell_ = std::vector<int>(1, 0);
    function_to_shell_ = std::vector<int>(1, 0);
    function_center_ = std::vector<int>(1, 0);
    shell_center_ = std::vector<int>(1, 0);
    center_to_nshell_ = std::vector<int>(1, 1);
    center_to_shell_ = std::vector<int>(1, 0);
    xyz_ = std::vector<double>(3, 0.0);
    puream_ = false;
    max_am_ = 0;
    max_nprimitive_ = 1;
    name_ = "(Empty Basis Set)";
    key_ = "(Empty Basis Set)";
    target_ = "(Empty Basis Set)";
    shells_[0] = GaussianShell(Gaussian, 0, nprimitive_, uoriginal_coefficients_.data(), ucoefficients_.data(),
                               uerd_coefficients_.data(), uexponents_.data(), GaussianType(0), 0, xyz_.data(), 0);
}

BasisSet::~BasisSet() {}

std::shared_ptr<BasisSet> BasisSet::build(std::shared_ptr<Molecule> /*molecule*/,
                                          const std::vector<ShellInfo> & /*shells*/) {
    // TRIAL//    //TODO fixme!!!
    // TRIAL//    auto basis = std::make_shared<BasisSet>();
    // TRIAL//    //    basis->molecule_ = molecule;
    // TRIAL//    //    basis->shells_ = shells;
    // TRIAL//    //    basis->refresh();
    // TRIAL//
    throw NOT_IMPLEMENTED_EXCEPTION();
    // TRIAL//    return basis;
}

void BasisSet::initialize_singletons() {
    if (initialized_shared_ == true) return;

    // Populate the exp_ao arrays
    for (int l = 0; l < LIBINT_MAX_AM; ++l) {
        for (int i = 0; i <= l; ++i) {
            int x = l - i;
            for (int j = 0; j <= i; ++j) {
                int y = i - j;
                int z = j;

                Vector3 xyz_ao(x, y, z);
                BasisSet::exp_ao[l].push_back(xyz_ao);
            }
        }
    }

    initialized_shared_ = true;
}

std::shared_ptr<Molecule> BasisSet::molecule() const { return molecule_; }

int BasisSet::n_ecp_core() const {
    int ncoreelectrons = 0;
    for (int A = 0; A < molecule_->natom(); A++) ncoreelectrons += n_ecp_core(molecule_->label(A));
    return ncoreelectrons;
}

static const std::vector<int> full_shell_values = {0, 2, 10, 18, 36, 54, 86, 118};

int BasisSet::atom_to_period(int Z) {
    if (Z > 118) {
        throw PSIEXCEPTION("Atomic number beyond Oganesson");
    }
    auto period = std::lower_bound(full_shell_values.begin(), full_shell_values.end(), Z) - full_shell_values.begin();
    return period;
}

int BasisSet::period_to_full_shell(int p) {
    if (p > 7) {
        throw PSIEXCEPTION("Atomic number beyond Oganesson");
    }
    return full_shell_values[p];
}

int BasisSet::n_frozen_core(const std::string &depth, SharedMolecule mol) {
    std::string local = depth;
    if (depth.empty()) local = Process::environment.options.get_str("FREEZE_CORE");

    SharedMolecule mymol = mol ? mol : molecule_;

    if (local == "FALSE" or local == "0") {
        return 0;
    } else if (local == "TRUE" or local == "1") {
        int num_frozen_el = 0;
        int mol_valence = -1 * mymol->molecular_charge();
        int largest_shell = 0;
        // Freeze the number of core electrons corresponding to the
        // nearest previous noble gas atom.  This means that the 4p block
        // will still have 3d electrons active.  Alkali earth atoms will
        // have one valence electron in this scheme.
        for (int A = 0; A < mymol->natom(); A++) {
            double Z = mymol->Z(A);
            // Exclude ghosted atoms from core-freezing
            if (Z > 0) {
                // Add ECPs to Z, the number of electrons less ECP-treated electrons
                int ECP = n_ecp_core(mymol->label(A));
                int current_shell = atom_to_period(Z + ECP);
                int delta = period_to_full_shell(current_shell - 1);
                // Keep track of the largest frozen shell, in case its a cationic species
                if (largest_shell < current_shell) {
                    largest_shell = current_shell;
                }
                // If center is a post-lanthanide or a post-actinide in Nth period,
                // freeze all 14 of its (N-2)f electrons too
                if (current_shell > 5) {
                    if ((Z + ECP - delta) >= 18) delta += 14;
                }
                // If this center has an ECP, some electrons are already frozen
                if (ECP > 0) delta -= ECP;
                // Keep track of current valence electrons
                mol_valence = mol_valence + Z - delta;
                num_frozen_el += delta;
            }
        }
        // If we are about to end up with no valence electrons,
        // unfreeze electrons from the largest shell in the molecule
        if (mol_valence <= 0)
            num_frozen_el -= period_to_full_shell(largest_shell - 1) - period_to_full_shell(largest_shell - 2);
        return num_frozen_el / 2;
    } else if (local == "POLICY") {
        // Look up what policy we've set
        std::vector<int> freeze_core_policy;
        freeze_core_policy = Process::environment.options.get_int_vector("FREEZE_CORE_POLICY");
        int nfzc = 0;
        for (int A = 0; A < mymol->natom(); A++) {
            // Exclude ghosted atoms from core-freezing
            double Z = mymol->Z(A);
            if (Z > 0) {
                // Add ECPs to Z, the number of electrons less ECP-treated electrons
                int true_Z = n_ecp_core(mymol->label(A)) + Z - 1;
                if (true_Z >= freeze_core_policy.size()) {
                    throw PSIEXCEPTION("Atomic number encountered greater than length of FREEZE_CORE_POLICY");
                }
                nfzc += freeze_core_policy[true_Z];
            }
        }
        return nfzc;
    } else {
        // Options are filtered in read_options.cc; allowed strings are:
        // TRUE, FALSE, -1, -2, -3
        int req_shell = -std::stoi(local, nullptr, 10);
        int num_frozen_el = 0;
        int mol_valence = -1 * mymol->molecular_charge();
        // Freeze the number of core electrons strictly corresponding to the
        // requested previous n-th shell.
        for (int A = 0; A < mymol->natom(); A++) {
            double Z = mymol->Z(A);
            // Exclude ghosted atoms from core-freezing
            if (Z > 0) {
                // Add ECPs to Z, the number of electrons less ECP-treated electrons
                int ECP = n_ecp_core(mymol->label(A));
                int current_shell = atom_to_period(Z + ECP);
                int delta = period_to_full_shell(std::max(current_shell - req_shell, 0));
                // If this center has an ECP, some electrons are already frozen
                if (delta < ECP)
                    throw PSIEXCEPTION(
                        "ECP on atom freezes more electrons than requested by choosing a previous shell.");
                if (ECP > 0) delta -= ECP;
                // Keep track of current valence electrons
                mol_valence = mol_valence + Z - delta;
                num_frozen_el += delta;
            }
        }
        // If we are about to end up with no valence electrons,
        // throw an exception.
        if (mol_valence <= 0) throw PSIEXCEPTION("Cannot freeze the requested previous shell: valence <= 0.");
        return num_frozen_el / 2;
    }
}

void BasisSet::print(std::string out) const {
    std::shared_ptr<psi::PsiOutStream> printer = (out == "outfile" ? outfile : std::make_shared<PsiOutStream>(out));
    printer->Printf("  Basis Set: %s\n", name_.c_str());
    printer->Printf("    Blend: %s\n", target_.c_str());
    printer->Printf("    Number of shells: %d\n", nshell());
    printer->Printf("    Number of basis functions: %d\n", nbf());
    printer->Printf("    Number of Cartesian functions: %d\n", nao());
    printer->Printf("    Spherical Harmonics?: %s\n", has_puream() ? "true" : "false");
    printer->Printf("    Max angular momentum: %d\n\n", max_am());
    if (has_ECP()) {
        printer->Printf("  Core potential: %s\n", name_.c_str());
        printer->Printf("    Number of shells: %d\n", n_ecp_shell());
        printer->Printf("    Number of ECP primitives: %d\n", n_ecp_primitive());
        printer->Printf("    Number of ECP core electrons: %d\n", n_ecp_core());
        printer->Printf("    Max angular momentum: %d\n\n", max_ecp_am());
    }
}

void BasisSet::print_by_level(std::string out, int level) const {
    if (level < 1)
        return;
    else if (level == 1)
        print(out);
    else if (level == 2)
        print_summary(out);
    else if (level > 2)
        print_detail(out);
}

void BasisSet::print_summary(std::string out) const {
    std::shared_ptr<psi::PsiOutStream> printer = (out == "outfile" ? outfile : std::make_shared<PsiOutStream>(out));

    printer->Printf("  -AO BASIS SET INFORMATION:\n");
    printer->Printf("    Name                   = %s\n", name_.c_str());
    printer->Printf("    Blend                  = %s\n", target_.c_str());
    printer->Printf("    Total number of shells = %d\n", nshell());
    printer->Printf("    Number of primitives   = %d\n", nprimitive_);
    printer->Printf("    Number of AO           = %d\n", nao_);
    printer->Printf("    Number of SO           = %d\n", nbf_);
    printer->Printf("    Maximum AM             = %d\n", max_am_);
    printer->Printf("    Spherical Harmonics    = %s\n", (puream_ ? "TRUE" : "FALSE"));
    printer->Printf("\n");

    printer->Printf("  -Contraction Scheme:\n");
    printer->Printf("    Atom   Type   All Primitives // Shells:\n");
    printer->Printf("   ------ ------ --------------------------\n");

    auto nprims = std::vector<int>(max_am_ + 1);
    auto nunique = std::vector<int>(max_am_ + 1);
    auto nshells = std::vector<int>(max_am_ + 1);
    auto amtypes = std::vector<char>(max_am_ + 1);

    for (int A = 0; A < molecule_->natom(); A++) {
        std::fill(nprims.begin(), nprims.end(), 0);
        std::fill(nunique.begin(), nunique.end(), 0);
        std::fill(nshells.begin(), nshells.end(), 0);

        printer->Printf("    %4d    ", A + 1);
        printer->Printf("%2s     ", molecule_->symbol(A).c_str());

        int first_shell = center_to_shell_[A];
        int n_shell = center_to_nshell_[A];

        for (int Q = 0; Q < n_shell; Q++) {
            const GaussianShell &shell = shells_[Q + first_shell];
            nshells[shell.am()]++;
            nunique[shell.am()] += shell.nprimitive();
            nprims[shell.am()] += shell.nprimitive();
            amtypes[shell.am()] = shell.amchar();
        }

        // All Primitives
        for (int l = 0; l < max_am_ + 1; l++) {
            if (nprims[l] == 0) continue;
            printer->Printf("%d%c ", nprims[l], amtypes[l]);
        }
        // Shells
        printer->Printf("// ");
        for (int l = 0; l < max_am_ + 1; l++) {
            if (nshells[l] == 0) continue;
            printer->Printf("%d%c ", nshells[l], amtypes[l]);
        }
        printer->Printf("\n");
    }
    printer->Printf("\n");

    if (has_ECP()) {
        printer->Printf("  -CORE POTENTIAL INFORMATION:\n");
        printer->Printf("    Total number of shells   = %d\n", n_ecp_shell());
        printer->Printf("    Number of terms          = %d\n", n_ecp_primitive_);
        printer->Printf("    Number of core electrons = %d\n", n_ecp_core());
        printer->Printf("    Maximum AM               = %d\n", max_ecp_am_);
        printer->Printf("\n");

        printer->Printf("  -Contraction Scheme:\n");
        printer->Printf("    Atom   Type  #elec    All terms // Shells:   \n");
        printer->Printf("   ------ ------ ----- --------------------------\n");

        auto nprims = std::vector<int>(max_ecp_am_ + 1);
        auto nunique = std::vector<int>(max_ecp_am_ + 1);
        auto nshells = std::vector<int>(max_ecp_am_ + 1);
        auto amtypes = std::vector<char>(max_ecp_am_ + 1);

        for (int A = 0; A < molecule_->natom(); A++) {
            memset((void *)nprims.data(), '\0', (max_ecp_am_ + 1) * sizeof(int));
            memset((void *)nunique.data(), '\0', (max_ecp_am_ + 1) * sizeof(int));
            memset((void *)nshells.data(), '\0', (max_ecp_am_ + 1) * sizeof(int));

            printer->Printf("    %4d    ", A + 1);
            printer->Printf("%2s     ", molecule_->symbol(A).c_str());

            int first_shell = center_to_ecp_shell_[A];
            int n_shell = center_to_ecp_nshell_[A];
            int ncoreelectrons = n_ecp_core(molecule_->label(A));
            printer->Printf("%2d   ", ncoreelectrons);

            if (ncoreelectrons) {
                for (int Q = 0; Q < n_shell; Q++) {
                    const GaussianShell &shell = ecp_shells_[Q + first_shell];
                    nshells[shell.am()]++;
                    nunique[shell.am()] += shell.nprimitive();
                    nprims[shell.am()] += shell.nprimitive();
                    amtypes[shell.am()] = shell.amchar();
                }

                // All Primitives
                for (int l = 0; l < max_ecp_am_ + 1; l++) {
                    if (nprims[l] == 0) continue;
                    printer->Printf("%d%c ", nprims[l], amtypes[l]);
                }
                // Shells
                if (n_shell) printer->Printf("// ");
                for (int l = 0; l < max_ecp_am_ + 1; l++) {
                    if (nshells[l] == 0) continue;
                    printer->Printf("%d%c ", nshells[l], amtypes[l]);
                }
            }
            printer->Printf("\n");
        }
        printer->Printf("\n");
    }
}

void BasisSet::print_detail(std::string out) const {
    print_summary(out);
    std::shared_ptr<psi::PsiOutStream> printer = (out == "outfile" ? outfile : std::make_shared<PsiOutStream>(out));

    printer->Printf("  ==> AO Basis Functions <==\n");
    printer->Printf("\n");
    printer->Printf("    [ %s ]\n", name_.c_str());
    if (has_puream())
        printer->Printf("    spherical\n");
    else
        printer->Printf("    cartesian\n");
    printer->Printf("    ****\n");

    for (int uA = 0; uA < molecule_->nunique(); uA++) {
        const int A = molecule_->unique(uA);

        printer->Printf("   %2s %3d\n", molecule_->symbol(A).c_str(), A + 1);

        int first_shell = center_to_shell_[A];
        int n_shell = center_to_nshell_[A];

        for (int Q = 0; Q < n_shell; Q++) shells_[Q + first_shell].print(out);

        printer->Printf("    ****\n");
    }

    printer->Printf("\n");
    if (has_ECP()) {
        printer->Printf("  ==> Core Potential Functions <==\n");
        printer->Printf("\n");
        printer->Printf("    [ %s ]\n", name_.c_str());
        printer->Printf("    ****\n");

        for (int uA = 0; uA < molecule_->nunique(); uA++) {
            const int A = molecule_->unique(uA);
            if (n_ecp_core(molecule_->label(A))) {
                int first_shell = center_to_ecp_shell_[A];
                int n_shell = center_to_ecp_nshell_[A];
                int shellam = 0;
                for (int shell = 0; shell < n_shell; ++shell)
                    shellam = ecp_shells_[shell + first_shell].am() > shellam ? ecp_shells_[shell + first_shell].am()
                                                                              : shellam;
                printer->Printf("   %2s %3d\n", molecule_->symbol(A).c_str(), A + 1);
                printer->Printf("   %2s-ECP  %d %3d\n", molecule_->symbol(A).c_str(), shellam,
                                n_ecp_core(molecule_->label(A)));

                for (int Q = 0; Q < n_shell; Q++) ecp_shells_[Q + first_shell].print(out);

                printer->Printf("    ****\n");
            }
        }

        printer->Printf("\n");
    }
}

std::string BasisSet::print_detail_cfour() const {
    char buffer[120];
    std::stringstream ss;
    std::string nameUpperCase = name_;
    to_upper(nameUpperCase);

    for (int uA = 0; uA < molecule_->nunique(); uA++) {
        const int A = molecule_->unique(uA);

        sprintf(buffer, "%s:P4_%d\n", molecule_->symbol(A).c_str(), A + 1);
        ss << buffer;
        sprintf(buffer, "Psi4 basis %s for element %s atom %d\n\n", nameUpperCase.c_str(), molecule_->symbol(A).c_str(),
                A + 1);
        ss << buffer;

        int first_shell = center_to_shell_[A];
        int n_shell = center_to_nshell_[A];

        int max_am_center = 0;
        for (int Q = 0; Q < n_shell; Q++)
            max_am_center =
                (shells_[Q + first_shell].am() > max_am_center) ? shells_[Q + first_shell].am() : max_am_center;

        std::vector<std::vector<int>> shell_per_am(max_am_center + 1);
        for (int Q = 0; Q < n_shell; Q++) shell_per_am[shells_[Q + first_shell].am()].push_back(Q);

        // Write number of shells in the basis set
        sprintf(buffer, "%3d\n", max_am_center + 1);
        ss << buffer;

        // Write angular momentum for each shell
        for (int am = 0; am <= max_am_center; am++) {
            sprintf(buffer, "%5d", am);
            ss << buffer;
        }
        sprintf(buffer, "\n");
        ss << buffer;

        // Write number of contracted basis functions for each shell
        for (int am = 0; am <= max_am_center; am++) {
            sprintf(buffer, "%5lu", shell_per_am[am].size());
            ss << buffer;
        }
        sprintf(buffer, "\n");
        ss << buffer;

        std::vector<std::vector<double>> exp_per_am(max_am_center + 1);
        std::vector<std::vector<double>> coef_per_am(max_am_center + 1);
        for (int am = 0; am <= max_am_center; am++) {
            // TODO: std::find safe on floats? seems to work
            // Collect unique exponents among all functions
            for (size_t Q = 0; Q < shell_per_am[am].size(); Q++) {
                for (int K = 0; K < shells_[shell_per_am[am][Q] + first_shell].nprimitive(); K++) {
                    if (!(std::find(exp_per_am[am].begin(), exp_per_am[am].end(),
                                    shells_[shell_per_am[am][Q] + first_shell].exp(K)) != exp_per_am[am].end())) {
                        exp_per_am[am].push_back(shells_[shell_per_am[am][Q] + first_shell].exp(K));
                    }
                }
            }

            // Collect coefficients for each exp among all functions, zero otherwise
            for (size_t Q = 0; Q < shell_per_am[am].size(); Q++) {
                for (size_t ep = 0, K = 0; ep < exp_per_am[am].size(); ep++) {
                    if (std::abs(exp_per_am[am][ep] - shells_[shell_per_am[am][Q] + first_shell].exp(K)) < 1.0e-8) {
                        coef_per_am[am].push_back(shells_[shell_per_am[am][Q] + first_shell].original_coef(K));
                        if ((K + 1) != (size_t)(shells_[shell_per_am[am][Q] + first_shell].nprimitive())) {
                            K++;
                        }
                    } else {
                        coef_per_am[am].push_back(0.0);
                    }
                }
            }
        }

        // Write number of exponents for each shell
        for (int am = 0; am <= max_am_center; am++) {
            sprintf(buffer, "%5lu", exp_per_am[am].size());
            ss << buffer;
        }
        sprintf(buffer, "\n\n");
        ss << buffer;

        for (int am = 0; am <= max_am_center; am++) {
            // Write exponents for each shell
            for (size_t ep = 0; ep < exp_per_am[am].size(); ep++) {
                if (exp_per_am[am][ep] >= 10000000.0)
                    sprintf(buffer, "%13.4f ", exp_per_am[am][ep]);
                else if (exp_per_am[am][ep] >= 1000000.0)
                    sprintf(buffer, "%13.5f ", exp_per_am[am][ep]);
                else if (exp_per_am[am][ep] >= 100000.0)
                    sprintf(buffer, "%13.6f ", exp_per_am[am][ep]);
                else
                    sprintf(buffer, "%14.7f", exp_per_am[am][ep]);
                ss << buffer;
                if (((ep + 1) % 5 == 0) || ((ep + 1) == exp_per_am[am].size())) {
                    sprintf(buffer, "\n");
                    ss << buffer;
                }
            }
            sprintf(buffer, "\n");
            ss << buffer;

            // Write contraction coefficients for each shell
            for (size_t ep = 0; ep < exp_per_am[am].size(); ep++) {
                for (size_t bf = 0; bf < shell_per_am[am].size(); bf++) {
                    sprintf(buffer, "%10.7f ", coef_per_am[am][bf * exp_per_am[am].size() + ep]);
                    ss << buffer;
                }
                sprintf(buffer, "\n");
                ss << buffer;
            }
            sprintf(buffer, "\n");
            ss << buffer;
        }
    }
    return ss.str();
}

const GaussianShell &BasisSet::shell(int si) const {
    if (si < 0 || si > nshell()) {
        outfile->Printf("BasisSet::shell(si = %d), requested a shell out-of-bound.\n", si);
        outfile->Printf("     Max shell size: %d\n", nshell());
        outfile->Printf("     Name: %s\n", name().c_str());
        throw PSIEXCEPTION("BasisSet::shell: requested shell is out-of-bounds.");
    }
    return shells_[si];
}

const GaussianShell &BasisSet::ecp_shell(int si) const {
    if (si < 0 || si > n_ecp_shell()) {
        outfile->Printf("BasisSet::ecp_shell(si = %d), requested a shell out-of-bound.\n", si);
        outfile->Printf("     Max shell size: %d\n", n_ecp_shell());
        outfile->Printf("     Name: %s\n", name().c_str());
        throw PSIEXCEPTION("BasisSet::ecp_shell: requested shell is out-of-bounds.");
    }
    return ecp_shells_[si];
}

const libint2::Shell &BasisSet::l2_shell(int si) const {
    if (si < 0 || si > nshell()) {
        outfile->Printf("Libint2 BasisSet::shell(si = %d), requested a shell out-of-bound.\n", si);
        outfile->Printf("     Max shell size: %d\n", nshell());
        outfile->Printf("     Name: %s\n", name().c_str());
        throw PSIEXCEPTION("BasisSet::shell: requested shell is out-of-bounds.");
    }
    return l2_shells_[si];
}

const GaussianShell &BasisSet::shell(int center, int si) const { return shell(center_to_shell_[center] + si); }

std::shared_ptr<BasisSet> BasisSet::zero_ao_basis_set() {
    // In the new implementation, we simply call the default constructor
    auto new_basis = std::make_shared<BasisSet>();
    return new_basis;
}

BasisSet::BasisSet(const std::string &basistype, SharedMolecule mol,
                   std::map<std::string, std::map<std::string, std::vector<ShellInfo>>> &shell_map,
                   std::map<std::string, std::map<std::string, std::vector<ShellInfo>>> &ecp_shell_map)
    : name_(basistype), molecule_(mol) {
    // Singletons
    initialize_singletons();

    int natom = molecule_->natom();

    /// These will tell us where the primitives for [basis][symbol] start and end, in the compact array
    std::map<std::string, std::map<std::string, int>> primitive_start;
    std::map<std::string, std::map<std::string, int>> primitive_end;
    std::map<std::string, std::map<std::string, int>> ecp_primitive_start;
    std::map<std::string, std::map<std::string, int>> ecp_primitive_end;

    /*
     * First, loop over the unique primitives, and store them
     */
    std::vector<double> uexps;
    std::vector<double> ucoefs;
    std::vector<double> uoriginal_coefs;
    std::vector<double> uerd_coefs;
    n_uprimitive_ = 0;
    std::map<std::string, std::map<std::string, std::vector<ShellInfo>>>::iterator basis_iter;
    for (basis_iter = shell_map.begin(); basis_iter != shell_map.end(); ++basis_iter) {
        const std::string &basis = basis_iter->first;
        std::map<std::string, std::vector<ShellInfo>> &symbol_map = shell_map[basis];
        std::map<std::string, std::vector<ShellInfo>>::iterator symbol_iter;
        for (symbol_iter = symbol_map.begin(); symbol_iter != symbol_map.end(); ++symbol_iter) {
            const std::string &label = symbol_iter->first;       // symbol --> label
            std::vector<ShellInfo> &shells = symbol_map[label];  // symbol --> label
            primitive_start[basis][label] = n_uprimitive_;       // symbol --> label
            for (size_t i = 0; i < shells.size(); ++i) {
                const ShellInfo &shell = shells[i];
                for (int prim = 0; prim < shell.nprimitive(); ++prim) {
                    uexps.push_back(shell.exp(prim));
                    ucoefs.push_back(shell.coef(prim));
                    uoriginal_coefs.push_back(shell.original_coef(prim));
                    uerd_coefs.push_back(shell.erd_coef(prim));
                    n_uprimitive_++;
                }
            }
            primitive_end[basis][label] = n_uprimitive_;  // symbol --> label
        }
    }

    /*
     * Now, loop over the ECP primitives, and store them
     */
    std::vector<double> uecpexps;
    std::vector<double> uecpcoefs;
    std::vector<int> uns;
    n_ecp_uprimitive_ = 0;
    for (basis_iter = ecp_shell_map.begin(); basis_iter != ecp_shell_map.end(); ++basis_iter) {
        const std::string &basis = basis_iter->first;
        std::map<std::string, std::vector<ShellInfo>> &symbol_map = ecp_shell_map[basis];
        std::map<std::string, std::vector<ShellInfo>>::iterator symbol_iter;
        for (symbol_iter = symbol_map.begin(); symbol_iter != symbol_map.end(); ++symbol_iter) {
            const std::string &label = symbol_iter->first;          // symbol --> label
            std::vector<ShellInfo> &shells = symbol_map[label];     // symbol --> label
            ecp_primitive_start[basis][label] = n_ecp_uprimitive_;  // symbol --> label
            for (size_t i = 0; i < shells.size(); ++i) {
                const ShellInfo &shell = shells[i];
                for (int prim = 0; prim < shell.nprimitive(); ++prim) {
                    uecpexps.push_back(shell.exp(prim));
                    uns.push_back(shell.nval(prim));
                    uecpcoefs.push_back(shell.coef(prim));
                    n_ecp_uprimitive_++;
                }
            }
            ecp_primitive_end[basis][label] = n_ecp_uprimitive_;  // symbol --> label
        }
    }

    /*
     * Count basis functions, shells and primitives
     */
    n_shells_ = 0;
    nprimitive_ = 0;
    nao_ = 0;
    nbf_ = 0;
    n_ecp_shells_ = 0;
    n_ecp_primitive_ = 0;
    for (int n = 0; n < natom; ++n) {
        const std::shared_ptr<CoordEntry> &atom = molecule_->atom_entry(n);
        std::string basis = atom->basisset(basistype);
        std::string label = atom->label();                         // symbol --> label
        std::vector<ShellInfo> &shells = shell_map[basis][label];  // symbol --> label
        for (size_t i = 0; i < shells.size(); ++i) {
            const ShellInfo &shell = shells[i];
            int nprim = shell.nprimitive();
            nprimitive_ += nprim;
            n_shells_++;
            nao_ += shell.ncartesian();
            nbf_ += shell.nfunction();
        }
        // ECP information
        if (!ecp_shell_map.empty()) {
            std::vector<ShellInfo> &ecpshells = ecp_shell_map[basis][label];
            for (size_t i = 0; i < ecpshells.size(); ++i) {
                const ShellInfo &ecpshell = ecpshells[i];
                int nprim = ecpshell.nprimitive();
                n_ecp_primitive_ += nprim;
                n_ecp_shells_++;
            }
        }
    }

    /*
     * Allocate arrays
     */
    n_prim_per_shell_ = std::vector<int>(n_shells_);
    // The unique primitives
    uexponents_ = std::vector<double>(n_uprimitive_);
    ucoefficients_ = std::vector<double>(n_uprimitive_);
    uoriginal_coefficients_ = std::vector<double>(n_uprimitive_);
    uerd_coefficients_ = std::vector<double>(n_uprimitive_);
    uecpexponents_ = std::vector<double>(n_ecp_uprimitive_);
    uecpcoefficients_ = std::vector<double>(n_ecp_uprimitive_);
    uecpns_ = std::vector<int>(n_ecp_uprimitive_);

    for (int i = 0; i < n_uprimitive_; ++i) {
        uexponents_[i] = uexps[i];
        ucoefficients_[i] = ucoefs[i];
        uoriginal_coefficients_[i] = uoriginal_coefs[i];
        uerd_coefficients_[i] = uerd_coefs[i];
    }
    for (int i = 0; i < n_ecp_uprimitive_; ++i) {
        uecpexponents_[i] = uecpexps[i];
        uecpcoefficients_[i] = uecpcoefs[i];
        uecpns_[i] = uns[i];
    }

    shell_first_ao_ = std::vector<int>(n_shells_);
    shell_first_exponent_ = std::vector<int>(n_shells_);
    shell_first_basis_function_ = std::vector<int>(n_shells_);
    shells_ = std::vector<GaussianShell>(n_shells_);
    ecp_shells_ = std::vector<GaussianShell>(n_ecp_shells_);
    ao_to_shell_ = std::vector<int>(nao_);
    function_to_shell_ = std::vector<int>(nbf_);
    function_center_ = std::vector<int>(nbf_);
    shell_center_ = std::vector<int>(n_shells_);
    center_to_nshell_ = std::vector<int>(natom);
    center_to_shell_ = std::vector<int>(natom);
    center_to_ecp_nshell_ = std::vector<int>(natom);
    center_to_ecp_shell_ = std::vector<int>(natom);
    ecp_shell_center_ = std::vector<int>(n_ecp_shells_);
    xyz_ = std::vector<double>(3 * natom);

    /*
     * Now loop over all atoms, and point to the appropriate unique data
     */
    int shell_count = 0;
    int ao_count = 0;
    int bf_count = 0;
    puream_ = false;
    max_am_ = 0;
    max_nprimitive_ = 0;
    for (int n = 0; n < natom; ++n) {
        const std::shared_ptr<CoordEntry> &atom = molecule_->atom_entry(n);
        std::string basis = atom->basisset(basistype);
        std::string label = atom->label();                         // symbol --> label
        std::vector<ShellInfo> &shells = shell_map[basis][label];  // symbol --> label
        int ustart = primitive_start[basis][label];                // symbol --> label
        int uend = primitive_end[basis][label];                    // symbol --> label
        int nshells = shells.size();
        center_to_nshell_[n] = nshells;
        center_to_shell_[n] = shell_count;
        Vector3 xyz = molecule_->xyz(n);
        xyz_[3 * n + 0] = xyz[0];
        xyz_[3 * n + 1] = xyz[1];
        xyz_[3 * n + 2] = xyz[2];
        int atom_nprim = 0;
        for (int i = 0; i < nshells; ++i) {
            const ShellInfo &thisshell = shells[i];
            ShellType shelltype = thisshell.shell_type();
            shell_first_ao_[shell_count] = ao_count;
            shell_first_basis_function_[shell_count] = bf_count;
            int shell_nprim = thisshell.nprimitive();
            int am = thisshell.am();
            max_nprimitive_ = shell_nprim > max_nprimitive_ ? shell_nprim : max_nprimitive_;
            max_am_ = max_am_ > std::abs(am) ? max_am_ : std::abs(am);
            shell_center_[shell_count] = n;
            GaussianType puream = thisshell.is_pure() ? Pure : Cartesian;
            if (puream) puream_ = true;
            if (shelltype == Gaussian) {
                // This is a regular Gaussian basis set
                shells_[shell_count] =
                    GaussianShell(shelltype, am, shell_nprim, &uoriginal_coefficients_[ustart + atom_nprim],
                                  &ucoefficients_[ustart + atom_nprim], &uerd_coefficients_[ustart + atom_nprim],
                                  &uexponents_[ustart + atom_nprim], puream, n, &xyz_.data()[3 * n], bf_count);
                shell_first_exponent_[shell_count] = ustart + atom_nprim;
                n_prim_per_shell_[shell_count] = shell_nprim;
            } else {
                throw PSIEXCEPTION("Unexpected shell type in BasisSet constructor!");
            }
            for (int thisbf = 0; thisbf < thisshell.nfunction(); ++thisbf) {
                function_to_shell_[bf_count] = shell_count;
                function_center_[bf_count++] = n;
            }
            for (int thisao = 0; thisao < thisshell.ncartesian(); ++thisao) {
                ao_to_shell_[ao_count++] = shell_count;
            }
            atom_nprim += shell_nprim;
            shell_count++;
        }
        if (atom_nprim != uend - ustart) {
            throw PSIEXCEPTION("Problem with nprimitive in basis set construction!");
        }
    }
    // Update the libint2 shell data
    update_l2_shells();

    /*
     * Now loop over ECPs and finalize metadata
     */
    max_ecp_am_ = -1;
    if (!ecp_shell_map.empty()) {
        int ecp_shell_count = 0;
        for (int n = 0; n < natom; ++n) {
            const std::shared_ptr<CoordEntry> &atom = molecule_->atom_entry(n);
            std::string basis = atom->basisset(basistype);
            std::string label = atom->label();                                 // symbol --> label
            std::vector<ShellInfo> &ecp_shells = ecp_shell_map[basis][label];  // symbol --> label
            int ustart = ecp_primitive_start[basis][label];                    // symbol --> label
            int uend = ecp_primitive_end[basis][label];                        // symbol --> label
            int n_ecp_shells = ecp_shells.size();
            center_to_ecp_nshell_[n] = n_ecp_shells;
            center_to_ecp_shell_[n] = ecp_shell_count;
            int atom_nprim = 0;
            for (int i = 0; i < n_ecp_shells; ++i) {
                const ShellInfo &thisshell = ecp_shells[i];
                ShellType shelltype = thisshell.shell_type();
                int ecp_shell_nprim = thisshell.nprimitive();
                int am = thisshell.am();
                max_ecp_am_ = max_ecp_am_ > std::abs(am) ? max_ecp_am_ : std::abs(am);
                ecp_shell_center_[ecp_shell_count] = n;
                if (shelltype == ECPType1 || shelltype == ECPType2) {
                    ecp_shells_[ecp_shell_count] = GaussianShell(
                        shelltype, am, ecp_shell_nprim, &uecpcoefficients_[ustart + atom_nprim],
                        &uecpexponents_[ustart + atom_nprim], &uecpns_[ustart + atom_nprim], n, &xyz_.data()[3 * n]);
                } else {
                    throw PSIEXCEPTION("Unknown ECP shell type in BasisSet constructor!");
                }
                atom_nprim += ecp_shell_nprim;
                ecp_shell_count++;
            }
            if (atom_nprim != uend - ustart) {
                throw PSIEXCEPTION("Problem with nprimitive in ECP basis set construction!");
            }
        }
    }
}

void BasisSet::update_l2_shells() {
    l2_shells_.resize(n_shells_);
    for (auto ishell = 0; ishell < n_shells_; ishell++) {
        auto am = shells_[ishell].am();
        auto center_index = shell_to_center(ishell);
        Vector3 xyz = molecule_->xyz(center_index);

        auto offset = shell_first_exponent_[ishell];
        auto nprim = n_prim_per_shell_[ishell];
        auto l2c = libint2::svector<double>(&uoriginal_coefficients_[offset], &uoriginal_coefficients_[offset + nprim]);
        auto l2e = libint2::svector<double>(&uexponents_[offset], &uexponents_[offset + nprim]);
        l2_shells_[ishell] = libint2::Shell{l2e, {{am, puream_, l2c}}, {{xyz[0], xyz[1], xyz[2]}}};
    }
}

std::string BasisSet::make_filename(const std::string &name) {
    // Modify the name of the basis set to generate a filename: STO-3G -> sto-3g
    std::string basisname = name;

    // First make it lower case
    std::transform(basisname.begin(), basisname.end(), basisname.begin(), ::tolower);

#if 0
    std::string format_underscore("_"); // empty string
    // Replace all '(' with '_'
    xpressive::sregex match_format = xpressive::as_xpr("(");
    basisname = regex_replace(basisname, match_format, format_underscore);

    // Replace all ')' with '_'
    match_format = xpressive::as_xpr(")");
    basisname = regex_replace(basisname, match_format, format_underscore);

    // Replace all ',' with '_'
    match_format = xpressive::as_xpr(",");
    basisname = regex_replace(basisname, match_format, format_underscore);

    // Replace all '*' with 's'
    match_format = xpressive::as_xpr("*");
    string format_star("s");
    basisname = regex_replace(basisname, match_format, format_star);

    // Replace all '+' with 'p'
    match_format = xpressive::as_xpr("+");
    string format_plus("p");
    basisname = regex_replace(basisname, match_format, format_plus);
#endif

    basisname = std::regex_replace(basisname, std::regex("\\(|\\)|,"), "_");
    basisname = std::regex_replace(basisname, std::regex("\\*"), "s");
    basisname = std::regex_replace(basisname, std::regex("\\+"), "p");

    // Add file extension
    basisname += ".gbs";

    return basisname;
}

void BasisSet::refresh() {
    // TODO FIXME!!!
    //    // Reset data to initial values
    //    nprimitive_ = 0;
    //    nao_ = 0;
    //    nbf_ = 0;
    //    max_am_ = 0;
    //    max_nprimitive_ = 0;
    //    puream_ = false;

    //    shell_first_basis_function_.clear(); shell_first_basis_function_.resize(nshell(), 0);
    //    shell_first_ao_.clear();             shell_first_ao_.resize(nshell(), 0);
    //    shell_center_.clear();               shell_center_.resize(nshell(), 0);
    //    function_center_.clear();
    //    center_to_nshell_.clear();           center_to_nshell_.resize(molecule_->natom(), 0);
    //    center_to_shell_.clear();            center_to_shell_.resize(molecule_->natom(), 0);
    //    center_to_shell_[0] = 0;

    //    int current_center = 0;

    //    for (int i=0; i<nshell(); ++i) {
    //        shell_center_[i]   = shells_[i].ncenter();
    //        shell_first_ao_[i] = nao_;
    //        shell_first_basis_function_[i] = nbf_;
    //        shells_[i].set_function_index(nbf_);

    //        center_to_nshell_[shell_center_[i]]++;
    //        if (current_center != shell_center_[i]) {
    //            center_to_shell_[shell_center_[i]] = i;
    //            current_center = shell_center_[i];
    //        }

    //        nprimitive_ += shells_[i].nprimitive();
    //        nao_        += shells_[i].ncartesian();
    //        nbf_        += shells_[i].nfunction();

    //        for (int m = 0; m < shells_[i].nfunction(); m++) {
    //            function_center_.push_back(shells_[i].ncenter());
    //        }

    //        if (max_am_ < shells_[i].am())
    //            max_am_ = shells_[i].am();

    //        if (max_nprimitive_ < shells_[i].nprimitive())
    //            max_nprimitive_ = shells_[i].nprimitive();

    //        if (puream_ == false && shells_[i].is_pure())
    //            puream_ = true;
    //    }

    //    function_to_shell_.resize(nbf());
    //    int ifunc = 0;
    //    for (int i=0; i<nshell(); ++i) {
    //        int nfun = shells_[i].nfunction();
    //        for (int j=0; j<nfun; ++j) {
    //            function_to_shell_[ifunc] = i;
    //            ifunc++;
    //        }
    //    }
    //    ao_to_shell_.resize(nao());
    //    ifunc = 0;
    //    for (int i=0; i<nshell(); ++i) {
    //        int nfun = shells_[i].ncartesian();
    //        for (int j=0; j<nfun; ++j) {
    //            ao_to_shell_[ifunc] = i;
    //            ifunc++;
    //        }
    //    }

    //    // Create a map that has a key/value pair
    //    // The key is the angular momentum function of the shell arranged in decending order
    //    // The value is the actual shell number
    //    typedef std::pair<int, int> am_to_shell_pair;
    //    std::multimap< int, int, std::less<int> > am_to_shell_list;
    //    for (int i=0; i < shells_.size(); i++) {
    //        am_to_shell_list.insert(am_to_shell_pair(shells_[i].nfunction(), i));
    //    }
    //    // This puts the sorted shell values into the sorted_shell_list_ vector
    //    // This can be used by the integral iterator to look up the value of the sorted shells
    //    std::multimap< int, int, std::less<int> >::iterator it;
    //    sorted_ao_shell_list_.clear();
    //    for (it=am_to_shell_list.begin(); it != am_to_shell_list.end(); it++) {
    //        //std::cout << "sorted shell size = " << it->first <<
    //        //        "\t, which belongs to shell number " << it->second << std::endl;
    //        sorted_ao_shell_list_.push_back(it->second);
    //    }
}

std::pair<std::vector<std::string>, std::shared_ptr<BasisSet>> BasisSet::test_basis_set(int /*max_am*/) {
    throw NOT_IMPLEMENTED_EXCEPTION();
#if 0
    int max_centers = 4;
    int max_primitives = 10;
    int max_shells;

    std::vector<int> nprim;
    nprim.push_back(10);
    nprim.push_back(1);
    nprim.push_back(6);
    nprim.push_back(1);
    nprim.push_back(2);
    nprim.push_back(1);
    nprim.push_back(1);
    nprim.push_back(1);
    nprim.push_back(1);
    nprim.push_back(1);

    std::vector<int> am;
    am.push_back(0);
    am.push_back(0);
    am.push_back(1);
    am.push_back(1);
    am.push_back(2);
    am.push_back(2);
    am.push_back(3);
    am.push_back(4);
    am.push_back(5);
    am.push_back(2);

    std::vector<std::vector<double> > c;
    c.push_back(std::vector<double>(10));
    c.push_back(std::vector<double>(1));
    c.push_back(std::vector<double>(6));
    c.push_back(std::vector<double>(1));
    c.push_back(std::vector<double>(2));
    c.push_back(std::vector<double>(1));
    c.push_back(std::vector<double>(1));
    c.push_back(std::vector<double>(1));
    c.push_back(std::vector<double>(1));
    c.push_back(std::vector<double>(1));

    std::vector<std::vector<double> > e;
    e.push_back(std::vector<double>(10));
    e.push_back(std::vector<double>(1));
    e.push_back(std::vector<double>(6));
    e.push_back(std::vector<double>(1));
    e.push_back(std::vector<double>(2));
    e.push_back(std::vector<double>(1));
    e.push_back(std::vector<double>(1));
    e.push_back(std::vector<double>(1));
    e.push_back(std::vector<double>(1));
    e.push_back(std::vector<double>(1));

    c[0][0] = 0.458878E-03;
    c[0][1] = 0.355070E-02;
    c[0][2] = 0.182618E-01;
    c[0][3] = 0.716650E-01;
    c[0][4] = 0.212346E+00;
    c[0][5] = 0.416203E+00;
    c[0][6] = 0.373020E+00;
    c[0][7] = 0.625054E-01;
    c[0][8] = 0.624532E-02;
    c[0][9] = 0.243374E-02;
    c[1][0] =     1.0;
    c[2][0] = 0.458878E-03;
    c[2][1] = 0.355070E-02;
    c[2][2] = 0.182618E-01;
    c[2][3] = 0.716650E-01;
    c[2][4] = 0.212346E+00;
    c[2][5] = 0.416203E+00;
    c[3][0] =     1.0;
    c[4][0] = 0.458878E-03;
    c[4][1] = 0.355070E-02;
    c[5][0] =     1.0;
    c[6][0] =     1.0;
    c[7][0] =     1.0;
    c[8][0] =     1.0;
    c[9][0] =     1.0;

    e[0][0] = 31700.0;
    e[0][1] =  4755.0;
    e[0][2] =  1082.0;
    e[0][3] =   306.0;
    e[0][4] =    99.0;
    e[0][5] =    33.0;
    e[0][6] =    13.0;
    e[0][7] =     4.0;
    e[0][8] =     2.0;
    e[0][9] =     0.5;
    e[1][0] =     1.0;
    e[2][0] = 31700.0;
    e[2][1] =  4755.0;
    e[2][2] =  1082.0;
    e[2][3] =   306.0;
    e[2][4] =    99.0;
    e[2][5] =    33.0;
    e[3][0] =     1.0;
    e[4][0] = 31700.0;
    e[4][1] =  4755.0;
    e[5][0] =     1.0;
    e[6][0] =     1.0;
    e[7][0] =     1.0;
    e[8][0] =     1.0;
    e[9][0] =     1.0;

    std::vector<std::string> labels;
    if (max_am > -1) {
        labels.push_back("S");
        labels.push_back("s");
        max_shells = 2;
    }
    labels.push_back("P");
    if (max_am > 0) {
        labels.push_back("p");
        max_shells = 4;
    }
    if (max_am > 1) {
        labels.push_back("D");
        labels.push_back("d");
        max_shells = 6;
    }
    if (max_am > 2) {
        labels.push_back("f");
        max_shells = 7;
    }
    if (max_am > 3) {
        labels.push_back("g");
        max_shells = 8;
    }
    if (max_am > 4) {
        labels.push_back("h");
        max_shells = 9;
    }
    if (max_am > 5) {
        labels.push_back("i");
        max_shells = 10;
    }

    auto new_basis = std::make_shared<BasisSet>();

    // Add 4 atoms to the molecule for this basis (max integal centers is 4 at the moment)
    new_basis->molecule_ = std::make_shared<Molecule>();
    // Ghost atoms are now handled differently, they are not added to the normal xyz information array,
    // but to the fxyz array.
    double x = 0.0;
    for (int A = 0; A < max_centers; A++) {
        new_basis->molecule_->add_atom(0, x, x, x);
        x += 1.0;
    }

    // Setup all the parameters needed for a zero basis set
    new_basis->shell_center_.resize(max_shells * max_centers);
    for (int A = 0; A < max_centers; A++)
        for (int Q = 0; Q < max_shells; Q++)
            new_basis->shell_center_[A*max_shells + Q] = A;
    new_basis->max_nprimitive_ = max_primitives;
    new_basis->max_am_ = max_am;

    // We'll time puream for now
    new_basis->puream_ = true;

    // Add shells
    for (int A = 0; A < max_centers; A++) {
        Vector3 center = new_basis->molecule_->fxyz(A);
        for (int Q = 0; Q < max_shells; Q++) {
            new_basis->shells_.push_back(GaussianShell(am[Q], c[Q], e[Q], Pure, A, center, 0));
        }
    }

    new_basis->refresh();

    return make_pair(labels, new_basis);
#endif
}

void BasisSet::move_atom(int atom, const Vector3 &trans) {
    int offset = 3 * atom;
    xyz_[offset + 0] += trans[0];
    xyz_[offset + 1] += trans[1];
    xyz_[offset + 2] += trans[2];
    for (int shell = 0; shell < nshell(); ++shell) {
        if (shells_[shell].ncenter() == atom) {
            l2_shells_[shell].O = std::array<double, 3>{xyz_[offset + 0], xyz_[offset + 1], xyz_[offset + 2]};
        }
    }
}

void BasisSet::compute_phi(double *phi_ao, double x, double y, double z) {
    zero_arr(phi_ao, nbf());

    int ao = 0;
    for (int ns = 0; ns < nshell(); ns++) {
        const GaussianShell &shell = shells_[ns];
        int am = shell.am();
        int nprim = shell.nprimitive();
        const double *a = shell.exps();
        const double *c = shell.coefs();

        const double *xyz = shell.center();
        double dx = x - xyz[0];
        double dy = y - xyz[1];
        double dz = z - xyz[2];
        double rr = dx * dx + dy * dy + dz * dz;

        double cexpr = 0;
        for (int np = 0; np < nprim; np++) cexpr += c[np] * exp(-a[np] * rr);

        if (puream_) {
            const auto s_transform = SphericalTransform(am);
            std::vector<double> cart_buffer(INT_NCART(am), 0.0);

            for (int l = 0; l < INT_NCART(am); l++) {
                Vector3 &components = exp_ao[am][l];
                cart_buffer[l] += pow(dx, static_cast<double>(components[0])) *
                                  pow(dy, static_cast<double>(components[1])) *
                                  pow(dz, static_cast<double>(components[2])) * cexpr;
            }

            for (int ind = 0; ind < s_transform.n(); ind++) {
                int lcart = s_transform.cartindex(ind);
                int lpure = s_transform.pureindex(ind);
                double coef = s_transform.coef(ind);

                phi_ao[ao + lpure] += coef * cart_buffer[lcart];
            }

        } else {
            for (int l = 0; l < INT_NCART(am); l++) {
                Vector3 &components = exp_ao[am][l];
                phi_ao[ao + l] += pow(dx, static_cast<double>(components[0])) *
                                  pow(dy, static_cast<double>(components[1])) *
                                  pow(dz, static_cast<double>(components[2])) * cexpr;
            }
        }

        ao += INT_NFUNC(puream_, am);
    }  // nshell
}
