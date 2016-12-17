/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2016 The Psi4 Developers.
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

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "psi4/dfep2/dfep2.h"
#include "psi4/libthce/thce.h"
#include "psi4/libthce/lreri.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libqt/qt.h"
#include "psi4/psifiles.h"
#include "psi4/libmints/vector.h"
#include "psi4/psi4-dec.h"
#include "psi4/libmints/matrix.h"

namespace psi { namespace dfep2 {

DFEP2Wavefunction::DFEP2Wavefunction(std::shared_ptr<Wavefunction> ref_wfn)
    : Wavefunction(Process::environment.options) {
    // Copy the wavefuntion then update
    shallow_copy(ref_wfn);

    conv_thresh_ = options_.get_double("EP2_CONVERGENCE");
    max_iter_ = options_.get_int("EP2_MAXITER");
    debug_ = options_.get_int("DEBUG");
    memory_doubles_ = (size_t) 0.07 * Process::environment.get_memory();

    dferi_ = DFERI::build(get_basisset("ORBITAL"), get_basisset("DF_BASIS_EP2"), options_);
    dferi_->print_header();

    AO_C_ = Ca_subset("AO", "ALL");
    AO_Cocc_ = Ca_subset("AO", "OCC");
    AO_Cvir_ = Ca_subset("AO", "VIR");
    AO_eps_ = epsilon_a_subset("AO", "ALL");

    for (size_t h = 0; h < nirrep_; h++) {
        for (size_t i = 0; i < nmopi_[h]; i++) {
            orbital_order_.push_back(std::tuple<double, size_t, size_t>(epsilon_a_->get(h, i), h, i));
        }
    }
    std::sort(orbital_order_.begin(), orbital_order_.end(),
              std::less<std::tuple<double, size_t, size_t>>());
}

std::vector<std::vector<std::pair<double, double>>> DFEP2Wavefunction::compute(std::vector<std::vector<size_t>> solve_orbs){

    // ==> Figure out which orbitals to compute <== /
    std::vector<std::tuple<size_t, size_t, size_t>> orb_positions;

    size_t nE = 0;
    size_t nfound = 0;
    for (size_t h = 0; h < solve_orbs.size(); h++){
        nE += solve_orbs[h].size();
    }

    if (solve_orbs.size() != nirrep_){
        throw PSIEXCEPTION("EP2: Size of solve_orbs does not match the number of irreps!");
    }

    for (size_t k = 0; k < orbital_order_.size(); k++){
        size_t h = std::get<1>(orbital_order_[k]);
        size_t current_i = std::get<2>(orbital_order_[k]);

        if (nfound == nE) break;

        for (const size_t& i: solve_orbs[h]){

            if (i == current_i){
                orb_positions.push_back(std::make_tuple(k, h, i));
                nfound++;
            }

            if (i > nmopi_[h]){
                throw PSIEXCEPTION("EP2: Orbital number is larger than the number of orbitals in the irrep!");
            }

        }

    }

    // Debug printing
    if (debug_ > 1) {
        printf("Nsolve %zu\n", nE);
        printf("Nfound %zu\n", nfound);
        for (size_t i = 0; i < orb_positions.size(); i++) {
            printf("orb_pos: %zu irrep: %zu irrep_pos: %zu\n", std::get<0>(orb_positions[i]),
                   std::get<1>(orb_positions[i]), std::get<2>(orb_positions[i]));
        }
    }

    if (nE != nfound){
        throw PSIEXCEPTION("EP2: AO -> SO mapping failed, did not find all orbitals to compute");
    }


    // ==> Set the epsilon vectors <== /
    size_t nocc = AO_Cocc_->colspi()[0];
    size_t nvir = AO_Cvir_->colspi()[0];

    std::vector<double> eps_E(nE);
    std::vector<double> denom_E(nE);
    std::vector<double> Enew(nE);
    double* AO_epsp = AO_eps_->pointer();
    for (size_t i = 0; i < orb_positions.size(); i++) {
        eps_E[i] = AO_epsp[std::get<0>(orb_positions[i])];
        denom_E[i] = eps_E[i];
        Enew[i] = eps_E[i];
    }

    std::vector<double> eps_occ(nocc);
    for (size_t i = 0; i < nocc; i++){
        eps_occ[i] = AO_epsp[i];
    }

    std::vector<double> eps_vir(nvir);
    for (size_t i = 0; i < nvir; i++){
        eps_vir[i] = AO_epsp[nocc + i];
    }


    // ==> Build the solve orbitals <== /
    SharedMatrix AO_CE(new Matrix("Solve orbitals", AO_C_->rowspi()[0], nE));

    double** AO_CEp = AO_CE->pointer();
    double** AO_Cp = AO_C_->pointer();
    int ce_start = 0;
    for (size_t k = 0; k < orb_positions.size(); k++) {
        int orb_copy = std::get<0>(orb_positions[k]);
        C_DCOPY(AO_C_->rowspi()[0], (AO_Cp[0] + orb_copy), AO_C_->colspi()[0],
                (AO_CEp[0] + ce_start), nE);
        ce_start++;
    }
    // AO_C_->print();
    // AO_Cocc_->print();
    // AO_Cvir_->print();
    // AO_CE->print();

    SharedMatrix C_Full = Matrix::horzcat({AO_Cocc_, AO_Cvir_, AO_CE});
    C_Full->set_name("Full C matrix");
    // C_Full->print();

    // ==> Transform DF integrals <== /

    dferi_->clear();
    dferi_->set_C(C_Full);
    dferi_->add_space("i", 0, nocc);
    dferi_->add_space("a", nocc, nocc + nvir);
    dferi_->add_space("E", nocc + nvir, nocc + nvir + nE);

    dferi_->add_pair_space("iaQ", "i", "a");
    dferi_->add_pair_space("aiQ", "i", "a", -1.0 / 2.0, true);
    dferi_->add_pair_space("EiQ", "E", "i");
    dferi_->add_pair_space("EaQ", "E", "a");
    dferi_->compute();

    std::map<std::string, std::shared_ptr<Tensor>>& dfints = dferi_->ints();
    size_t nQ = dferi_->size_Q();

    // ==> Setup contraction <== /

    std::vector<double> Eold(nE);
    std::vector<double> Esigma(nE);
    std::vector<double> Ederiv(nE);
    std::vector<double> Eerror(nE);

    size_t fstat;

    // Always read
    std::shared_ptr<Tensor> EaQT = dfints["EaQ"];
    SharedMatrix EaQ(new Matrix("EaQ", nE * nvir, nQ));
    double* EaQp = EaQ->pointer()[0];
    FILE* EaQF = EaQT->file_pointer();
    fseek(EaQF, 0L, SEEK_SET);
    fstat = fread(EaQp, sizeof(double), nE * nvir * nQ, EaQF);

    // Read blocks
    std::shared_ptr<Tensor> aiQT = dfints["aiQ"];
    SharedMatrix aiQ(new Matrix("aiQ", nocc * nvir, nQ));
    double* aiQp = aiQ->pointer()[0];
    FILE* aiQF = aiQT->file_pointer();
    fseek(aiQF, 0L, SEEK_SET);
    fstat = fread(aiQp, sizeof(double), nocc * nvir * nQ, aiQF);
    aiQ->print();

    SharedMatrix I_Evvo = Matrix::doublet(EaQ, aiQ, false, true);

    // Always read
    std::shared_ptr<Tensor> EiQT = dfints["EiQ"];
    SharedMatrix EiQ(new Matrix("EiQ", nE * nocc, nQ));
    double* EiQp = EiQ->pointer()[0];
    FILE* EiQF = EiQT->file_pointer();
    fseek(EiQF, 0L, SEEK_SET);
    fstat = fread(EiQp, sizeof(double), nE * nocc * nQ, EiQF);

    // Read blocks
    std::shared_ptr<Tensor> iaQT = dfints["iaQ"];
    SharedMatrix iaQ(new Matrix("iaQ", nocc * nvir, nQ));
    double* iaQp = iaQ->pointer()[0];
    FILE* iaQF = iaQT->file_pointer();
    fseek(iaQF, 0L, SEEK_SET);
    fstat = fread(iaQp, sizeof(double), nocc * nvir * nQ, iaQF);

    SharedMatrix I_Eoov = Matrix::doublet(EiQ, iaQ, false, true);

    // SharedMatrix tmp_mat(new Matrix(nE * nvir, nvir * nocc));
    // SharedMatrix eps_mat(new Matrix(nE * nvir, nvir * nocc));

    // ==> Iterate <== /

    for (size_t iter = 0; iter < max_iter_; iter++) {
        for (size_t i = 0; i < nE; i++) {
            Esigma[i] = 0.0;
            Ederiv[i] = 0.0;
            Eold[i] = Enew[i];
        }

        // => Excitations <= //
        // sigma <= (Eabi - Ebai) * Eabi / (E - v - v + o)
        double** I_Evvop = I_Evvo->pointer();
        // double** tmpp = tmp_mat->pointer();
        // double** epsp = eps_mat->pointer();
        for (size_t e = 0; e < nE; e++) {
            // printf("%16.8f\n", denom_E[e]);
            for (size_t a = 0; a < nvir; a++) {
                for (size_t b = 0; b < nvir; b++) {
                    for (size_t i = 0; i < nocc; i++) {
                        double Eabi = I_Evvop[e * nvir + a][b * nocc + i];
                        double Ebai = I_Evvop[e * nvir + b][a * nocc + i];
                        double numer = (2.0 * Eabi - Ebai) * Eabi;
                        double denom = (denom_E[e] - eps_vir[a] - eps_vir[b] + eps_occ[i]);
                        // tmpp[e * nvir + a][b * nocc + i] = numer;
                        // epsp[e * nvir + a][b * nocc + i] = denom;

                        Esigma[e] += numer / denom;
                        Ederiv[e] += numer / (denom * denom);
                    }
                }
            }
        }

        // => De-excitations <= //
        // sigma <= (Eija - Ejia) * Eija / (E - o - o + v)
        double** I_Eoovp = I_Eoov->pointer();
        for (size_t e = 0; e < nE; e++) {
            for (size_t i = 0; i < nocc; i++) {
                for (size_t j = 0; j < nocc; j++) {
                    for (size_t a = 0; a < nvir; a++) {
                        double Eija = I_Eoovp[e * nocc + i][j * nvir + a];
                        double Ejia = I_Eoovp[e * nocc + j][i * nvir + a];
                        double numer = (2.0 * Eija - Ejia) * Eija;
                        double denom = (denom_E[e] - eps_occ[i] - eps_occ[j] + eps_vir[a]);
                        // tmpp[e * nocc + i][j * nvir + a] = numer;
                        // epsp[e * nocc + j][i * nvir + a] = denom;

                        Esigma[e] += numer / denom;
                        Ederiv[e] += numer / (denom * denom);
                    }
                }
            }
        }

        double max_error = 0.0;
        double mean_error = 0.0;
        size_t nremain = nE;
        for (size_t i = 0; i < nE; i++) {
            // printf("%16.8f  ", Enew[i]);

            // Compute new energy and update
            Enew[i] = eps_E[i] + Esigma[i];
            denom_E[i] = Eold[i] - ((Eold[i] - Enew[i]) / (1 + Ederiv[i]));

            // Compute stats
            Eerror[i] = std::fabs(Enew[i] - Eold[i]);
            mean_error += Eerror[i];
            if (Eerror[i] < conv_thresh_) {
                nremain--;
            }
            if (Eerror[i] > max_error) {
                max_error = Eerror[i];
            }

            // Reset E's
            Eold[i] = Enew[i];
            Esigma[i] = 0.0;
            Ederiv[i] = 0.0;
        }
        mean_error /= (size_t)nE;
        // printf("\n");

        printf("%zu %16.8f %16.8f %zu\n", (iter + 1), mean_error, max_error, nremain);

        if (max_error < conv_thresh_) break;
    }

    // Build output array and remap symmetry
    std::vector<std::vector<std::pair<double, double>>> ret;
    for (size_t h = 0; h < nirrep_; h++) {
        std::vector<std::pair<double, double>> tmp(solve_orbs[h].size());
        ret.push_back(tmp);
    }

    for (size_t k = 0; k < orb_positions.size(); k++) {
        size_t h = std::get<1>(orb_positions[k]);
        size_t i = std::get<2>(orb_positions[k]);
        ret[h][i].first = Enew[k];
        ret[h][i].second = Eerror[k];
    }

    return ret;

}

}}
