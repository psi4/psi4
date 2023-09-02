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

/*
 * frac.cc.
 *
 * frac handles SCF with fractional occupation numbers. Each fractionally occupied orbital is regarded as occupied
 * for purposes of docc, socc, nalpha, nbeta, etc. However, each orbital's contribution to the density matrix is
 * scaled by the fractional occupation number. The density matrix no longer describes a single Slater determinant
 * but an ensemble of Slater determinants. This is useful for finite-temperature theory and removing artifacts
 * where some orbitals but not others of a degenerate set are held occupied, artificially breaking the degeneracy.
 *
 * The way this works:
 * 1) form_C returns 1's normalized C matrices.
 * 2) find_occupation determines the occupations of said matrix via either Aufbau or MOM selection
 * 3) find_occupation then calls this, which renormalizes the C matrices with \sqrt(val) for each frac occ
 *    These renormalized C matrices produce the desired density matrices, but are not orthonormal.
 * 4) find_D is then computed with the renormalized C matrices, and everything is transparent until 1) on the next cycle
 * 5) frac_renormalize sets the wfn's C matrices to be un-renormalized.
 *
 * Some executive decisions:
 *  -DIIS: Upon FRAC start, the old DIIS info is nuked, and DIIS_START is incremented by the current iteration count.
 *         Thus, DIIS begins again on the next iteration. The exception is if you set FRAC_DIIS to false,
 *         in which case DIIS will cease for good upon FRAC start.
 *  -MOM:  To use MOM with FRAC, set MOM_START to either FRAC_START or FRAC_START+1 (I think I prefer the latter).
 */

#include "hf.h"

#include "psi4/liboptions/liboptions.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/integral.h"
#include "psi4/libmints/sointegral_onebody.h"
#include "psi4/libmints/factory.h"
#include "psi4/libqt/qt.h"
#include "psi4/libpsi4util/PsiOutStream.h"

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <vector>
#include <utility>
#include <tuple>

namespace psi {
namespace scf {

void HF::frac() {
    // Perhaps no frac?
    if (iteration_ < options_.get_int("FRAC_START") || options_.get_int("FRAC_START") == 0) return;
    frac_performed_ = true;

    // First frac iteration, blow away the diis and print the frac task
    if (iteration_ == options_.get_int("FRAC_START")) {
        if (!(options_.get_str("REFERENCE") == "UHF" || options_.get_str("REFERENCE") == "UKS"))
            throw PSIEXCEPTION("Fractional Occupation SCF is only implemented for UHF/UKS");

        if (!options_["FRAC_OCC"].size())
            throw PSIEXCEPTION("Fractional Occupation SCF requested, but empty FRAC_OCC/FRAC_VAL vector");

        if (options_["FRAC_OCC"].size() != options_["FRAC_VAL"].size())
            throw PSIEXCEPTION("Fractional Occupation SCF: FRAC_OCC/FRAC_VAL are of different dimensions");

        if (input_docc_ || input_socc_) throw PSIEXCEPTION("Fractional Occupation SCF: Turn off DOCC/SOCC");

        if (options_.get_int("MOM_START") <= options_.get_int("FRAC_START") && options_.get_int("MOM_START") != 0)
            throw PSIEXCEPTION("Fractional Occupation SCF: MOM must start after FRAC");

        if (MOM_excited_) throw PSIEXCEPTION("Fractional Occupation SCF: Don't try an excited-state MOM");

        // Close off a previous burn-in SCF
        outfile->Printf("\n");
        print_orbitals();

        // frac header
        outfile->Printf("\n  ==> Fractionally-Occupied SCF Iterations <==\n\n");
        for (int ind = 0; ind < options_["FRAC_OCC"].size(); ind++) {
            int i = options_["FRAC_OCC"][ind].to_integer();
            double val = options_["FRAC_VAL"][ind].to_double();

            // Throw if user requests frac occ above nalpha/nbeta
            int max_i = (i > 0 ? nalpha_ : nbeta_);
            if (std::abs(i) > max_i) {
                if (i > 0)
                    nalpha_++;
                else
                    nbeta_++;
            }

            if (val < 0.0)
                throw PSIEXCEPTION("Fractional Occupation SCF: Occupations must be non-negative.");

            outfile->Printf("    %-5s orbital %4d will contain %11.3E electron.\n", (i > 0 ? "Alpha" : "Beta"),
                            std::abs(i), val);
        }
        outfile->Printf("\n");

        // Make sure diis restarts correctly/frac plays well with MOM
        if (initialized_diis_manager_) {
            diis_manager_.attr("delete_diis_file")();
            diis_manager_ = py::none();
            initialized_diis_manager_ = false;
            diis_start_ += iteration_ + 1;
        }

        // Turn yonder DIIS off if requested
        if (!options_.get_bool("FRAC_DIIS")) {
            diis_enabled_ = false;
        }

        // Load the old orbitals in if requested
        if (options_.get_bool("FRAC_LOAD")) {
            outfile->Printf("    Orbitals reloaded from file, your previous iterations are garbage.\n\n");
            throw PSIEXCEPTION("FRAC_LOAD is currently not an available feature");
            // load_orbitals();
        }

        // Prevent spurious convergence (technically this iteration comes from the N-electron system anyways)
        frac_performed_ = false;
    }

    // Every frac iteration: renormalize the Ca/Cb matrices
    frac_helper();
}

void HF::frac_renormalize() {
    if ((options_.get_int("FRAC_START") == 0)       // frac disabled
        || !options_.get_bool("FRAC_RENORMALIZE"))  // || don't renormalize C
        return;

    // Renormalize the fractional occupations back to 1, if possible before storage
    // We need to cache the old quantities because the renormalization isn't invertible
    // for an occupation of 0. A user may want this during, say, frac_traverse.
    outfile->Printf("    FRAC: Renormalizing orbitals to 1.0 for storage.\n\n");

    Ca_ = unscaled_Ca_;
    Cb_ = unscaled_Cb_;
}


void HF::frac_helper() {

    // Save the orbitals, without scaling factors
    unscaled_Ca_ = Ca_->clone();
    unscaled_Cb_ = Cb_->clone();

    // Sort the eigenvalues in the usual manner
    std::vector<std::tuple<double, int, int> > pairs_a;
    std::vector<std::tuple<double, int, int> > pairs_b;
    for (int h = 0; h < epsilon_a_->nirrep(); ++h) {
        for (int i = 0; i < epsilon_a_->dimpi()[h]; ++i)
            pairs_a.push_back(std::tuple<double, int, int>(epsilon_a_->get(h, i), h, i));
    }
    for (int h = 0; h < epsilon_b_->nirrep(); ++h) {
        for (int i = 0; i < epsilon_b_->dimpi()[h]; ++i)
            pairs_b.push_back(std::tuple<double, int, int>(epsilon_b_->get(h, i), h, i));
    }
    sort(pairs_a.begin(), pairs_a.end());
    sort(pairs_b.begin(), pairs_b.end());

    // Renormalize the C matrix entries
    for (int ind = 0; ind < options_["FRAC_OCC"].size(); ind++) {
        int i = options_["FRAC_OCC"][ind].to_integer();
        double val = options_["FRAC_VAL"][ind].to_double();
        bool is_alpha = (i > 0);
        i = std::abs(i) - 1;  // Back to C ordering

        int h; double** Cp; int j;
        if (is_alpha) {
            h = std::get<1>(pairs_a[i]);
            j = std::get<2>(pairs_a[i]);
            Cp = Ca_->pointer(h);
        } else {
            h = std::get<1>(pairs_b[i]);
            j = std::get<2>(pairs_b[i]);
            Cp = Cb_->pointer(h);
        }

        int nso = Ca_->rowspi()[h];
        int nmo = Ca_->colspi()[h];

        C_DSCAL(nso, std::sqrt(val), &Cp[0][j], nmo);
    }
}

void HF::compute_spin_contamination() {
    if (!(options_.get_str("REFERENCE") == "UHF" || options_.get_str("REFERENCE") == "UKS" ||
          options_.get_str("REFERENCE") == "CUHF"))
        return;

    auto nalpha = (double)nalpha_;
    auto nbeta = (double)nbeta_;

    // Adjust for fractional occupation
    if (frac_performed_) {
        for (int ind = 0; ind < options_["FRAC_OCC"].size(); ind++) {
            int i = options_["FRAC_OCC"][ind].to_integer();
            double val = options_["FRAC_VAL"][ind].to_double();
            bool is_alpha = (i > 0);
            if (is_alpha) {
                nalpha -= (1.0 - val);
            } else {
                nbeta -= (1.0 - val);
            }
        }
    }

    SharedMatrix S = SharedMatrix(factory_->create_matrix("S (Overlap)"));
    auto fact = std::make_shared<IntegralFactory>(basisset_, basisset_, basisset_, basisset_);
    std::shared_ptr<OneBodySOInt> so_overlap(fact->so_overlap());
    so_overlap->compute(S);

    double dN = 0.0;

    for (int h = 0; h < S->nirrep(); h++) {
        int nbf = S->colspi()[h];
        int nmo = Ca_->colspi()[h];
        int na = nalphapi_[h];
        int nb = nbetapi_[h];
        if (na == 0 || nb == 0 || nbf == 0 || nmo == 0) continue;

        auto Ht = std::make_shared<Matrix>("H Temp", nbf, nb);
        auto Ft = std::make_shared<Matrix>("F Temp", na, nb);

        double** Sp = S->pointer(h);
        double** Cap = Ca_->pointer(h);
        double** Cbp = Cb_->pointer(h);
        double** Htp = Ht->pointer(0);
        double** Ftp = Ft->pointer(0);

        C_DGEMM('N', 'N', nbf, nb, nbf, 1.0, Sp[0], nbf, Cbp[0], nmo, 0.0, Htp[0], nb);
        C_DGEMM('T', 'N', na, nb, nbf, 1.0, Cap[0], nmo, Htp[0], nb, 0.0, Ftp[0], nb);

        dN += C_DDOT(na * (long int)nb, Ftp[0], 1, Ftp[0], 1);
    }

    double nmin = (nbeta < nalpha ? nbeta : nalpha);
    double dS = nmin - dN;
    double nm = (nalpha - nbeta) / 2.0;
    double S2 = std::fabs(nm) * (std::fabs(nm) + 1.0);

    outfile->Printf("   @Spin Contamination Metric: %17.9E\n", dS);
    outfile->Printf("   @S^2 Expected:              %17.9E\n", S2);
    outfile->Printf("   @S^2 Observed:              %17.9E\n", S2 + dS);
    outfile->Printf("   @S   Expected:              %17.9E\n", nm);
    outfile->Printf("   @S   Observed:              %17.9E\n", nm);

    if (frac_performed_) {
        outfile->Printf("   @Nalpha:                    %17.9E\n", nalpha);
        outfile->Printf("   @Nbeta:                     %17.9E\n", nbeta);
    }
    outfile->Printf("\n");
}
}  // namespace scf
}  // namespace psi
