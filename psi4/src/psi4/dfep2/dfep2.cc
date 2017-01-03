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
#include <iomanip>
#include "psi4/dfep2/dfep2.h"
#include "psi4/libthce/thce.h"
#include "psi4/libthce/lreri.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libqt/qt.h"
#include "psi4/psifiles.h"
#include "psi4/libmints/vector.h"
#include "psi4/psi4-dec.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/libpsio/psio.h"
#include "psi4/libpsio/aiohandler.h"

#ifdef _OPENMP
#include <omp.h>
#endif

namespace psi { namespace dfep2 {

DFEP2Wavefunction::DFEP2Wavefunction(std::shared_ptr<Wavefunction> ref_wfn)
    : Wavefunction(Process::environment.options) {
    // Copy the wavefuntion then update
    shallow_copy(ref_wfn);

    outfile->Printf("\n");
    outfile->Printf("         ---------------------------------------------------------\n");
    outfile->Printf("                  Second-Order Electron Propagator Theory\n");
    outfile->Printf("\n");
    outfile->Printf("                            by Daniel G. A. Smith\n");
    outfile->Printf("         ---------------------------------------------------------\n\n");


    conv_thresh_ = options_.get_double("EP2_CONVERGENCE");
    max_iter_ = options_.get_int("EP2_MAXITER");
    debug_ = options_.get_int("DEBUG");
    memory_doubles_ = (size_t)(0.1 * (double)Process::environment.get_memory());

    dferi_ = DFERI::build(get_basisset("ORBITAL"), get_basisset("DF_BASIS_EP2"), options_);
    dferi_->print_header();

    AO_C_ = Ca_subset("AO", "ALL");
    AO_Cocc_ = Ca_subset("AO", "OCC");
    AO_Cvir_ = Ca_subset("AO", "VIR");
    AO_eps_ = epsilon_a_subset("AO", "ALL");

    unit_ = PSIF_DFMP2_QIA;

    for (size_t h = 0; h < nirrep_; h++) {
        for (size_t i = 0; i < nmopi_[h]; i++) {
            orbital_order_.push_back(std::tuple<double, size_t, size_t>(epsilon_a_->get(h, i), h, i));
        }
    }
    std::sort(orbital_order_.begin(), orbital_order_.end(),
              std::less<std::tuple<double, size_t, size_t>>());

    num_threads_ = 1;
    #ifdef _OPENMP
        num_threads_ = omp_get_max_threads();
    #endif
}

std::vector<std::vector<std::pair<double, double>>> DFEP2Wavefunction::compute(std::vector<std::vector<size_t>> solve_orbs){

    // ==> Figure out which orbitals to compute <== /
    std::vector<std::tuple<size_t, size_t, size_t, size_t>> orb_positions;

    size_t nE = 0;
    size_t nfound = 0;
    for (size_t h = 0; h < solve_orbs.size(); h++){
        nE += solve_orbs[h].size();
    }

    if (solve_orbs.size() != nirrep_){
        throw PSIEXCEPTION("EP2: Size of solve_orbs does not match the number of irreps!");
    }

    std::vector<size_t> h_rel_index(nirrep_);

    for (size_t k = 0; k < orbital_order_.size(); k++){
        size_t h = std::get<1>(orbital_order_[k]);
        size_t current_i = std::get<2>(orbital_order_[k]);

        if (nfound == nE) break;

        for (const size_t& i: solve_orbs[h]){

            if (i == current_i){

                orb_positions.push_back(std::make_tuple(k, h, i, h_rel_index[h]));
                h_rel_index[h]++;
                nfound++;
            }

            if (i > nmopi_[h]){
                throw PSIEXCEPTION("EP2: Orbital number is larger than the number of orbitals in the irrep!");
            }

        }

    }
    outfile->Printf("  ==> Algorithm <==\n\n");
    outfile->Printf("  Number of orbitals:     %8zu\n", nE);
    outfile->Printf("  Maximum iterations:     %8zu\n", max_iter_);
    outfile->Printf("  Convergence threshold:  %3.2e\n\n", conv_thresh_);

    // Debug printing
    if (debug_ > 1) {
        outfile->Printf("\n\nNsolve %zu\n", nE);
        outfile->Printf("Nfound %zu\n", nfound);
        for (size_t i = 0; i < orb_positions.size(); i++) {
            outfile->Printf("orb_pos: %zu irrep: %zu irrep_pos: %zu\n",
                            std::get<0>(orb_positions[i]), std::get<1>(orb_positions[i]),
                            std::get<2>(orb_positions[i]));
        }
        outfile->Printf("\n\n");
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
    // dferi_->add_pair_space("aiQ", "i", "a", -0.5, true);
    dferi_->add_pair_space("iEQ", "i", "E");
    dferi_->add_pair_space("aEQ", "a", "E");
    dferi_->compute();

    std::map<std::string, std::shared_ptr<Tensor>>& dfints = dferi_->ints();
    size_t nQ = dferi_->size_Q();

    // ==> Build ERI's <== /

    std::shared_ptr<PSIO> psio = PSIO::shared_object();
    std::shared_ptr<AIOHandler> aio(new AIOHandler(psio));

    psio->open(unit_, PSIO_OPEN_OLD);

    aio->zero_disk(unit_, "EP2 I_ovvE Integrals", (ULI)(nocc * nvir), (ULI)(nvir * nE));
    aio->zero_disk(unit_, "EP2 I_vooE Integrals", (ULI)(nocc * nvir), (ULI)(nocc * nE));
    aio->synchronize();

    // How much memory are we working with?
    size_t E_tensor_size = nE * nvir * nQ + nE * nocc * nQ;
    size_t I_block_sizes = nvir * nvir * nE + nocc * nvir * nE + nvir * nQ;
    size_t free_doubles = memory_doubles_ - E_tensor_size;

    size_t block_size = free_doubles / I_block_sizes;
    if (block_size > nocc) block_size = nocc;
    // block_size = 2;

    size_t nblocks = 1 + ((nocc - 1) / block_size);

    if (debug_ > 1){
        outfile->Printf("\n\n");
        outfile->Printf("memory_doubles_ %zu\n", memory_doubles_);
        outfile->Printf("free_doubles    %zu\n", free_doubles);
        outfile->Printf("E_tensor_size   %zu\n", E_tensor_size);
        outfile->Printf("I_block_sizes   %zu\n", I_block_sizes);
        outfile->Printf("Block size      %zu\n", block_size);
        outfile->Printf("N Block         %zu\n", nblocks);
        outfile->Printf("\n\n");
    }

    if (block_size < 1){
        std::stringstream message;
        double mem_gb = ((double)(E_tensor_size + I_block_sizes) / 1.e10);
        message << "DF-EP2 requires at least is nvir^2 * number solve orbitals in memory." << std::endl;
        message << "       After taxes this is " << std::setprecision(2) << mem_gb << " GB of memory.";

        throw PSIEXCEPTION(message.str());
    }

    size_t fstat;


    // Read in part of the tensors
    std::shared_ptr<Tensor> aEQT = dfints["aEQ"];
    SharedMatrix aEQ(new Matrix("aEQ", nE * nvir, nQ));
    double* aEQp = aEQ->pointer()[0];
    FILE* aEQF = aEQT->file_pointer();
    fseek(aEQF, 0L, SEEK_SET);
    fstat = fread(aEQp, sizeof(double), nE * nvir * nQ, aEQF);

    std::shared_ptr<Tensor> iEQT = dfints["iEQ"];
    SharedMatrix iEQ(new Matrix("iEQ", nE * nocc, nQ));
    double* iEQp = iEQ->pointer()[0];
    FILE* iEQF = iEQT->file_pointer();
    fseek(iEQF, 0L, SEEK_SET);
    fstat = fread(iEQp, sizeof(double), nE * nocc * nQ, iEQF);

    // Allocate temps
    SharedMatrix block_iaQ(new Matrix(block_size * nvir, nQ));
    SharedMatrix temp_ovvE(new Matrix(block_size * nvir, nvir * nE));
    SharedMatrix temp_ovoE(new Matrix(block_size * nvir, nocc * nE));

    // Pointer to read from
    std::shared_ptr<Tensor> iaQT = dfints["iaQ"];
    FILE* iaQF = iaQT->file_pointer();
    double* block_iaQp = block_iaQ->pointer()[0];
    fseek(iaQF, 0L, SEEK_SET);

    psio_address ovvE_addr = psio_get_address(PSIO_ZERO, 0);
    psio_address ovoE_addr = psio_get_address(PSIO_ZERO, 0);
    psio_address vooE_addr = psio_get_address(PSIO_ZERO, 0);

    for (size_t block = 0; block < nblocks; block++){
        size_t bstart = block * block_size;

        // Make sure we dont read past the line
        if ((bstart + block_size) > nocc){
            block_size = nocc - block * block_size;
            block_iaQ->zero();
        }

        // Read a IA block
        fstat = fread(block_iaQp, sizeof(double), block_size * nvir * nQ, iaQF);

        // Write out OVVE
        temp_ovvE->gemm(false, true, 1.0, block_iaQ, aEQ, 0.0);
        psio_->write(unit_, "EP2 I_ovvE Integrals", (char*)temp_ovvE->pointer()[0],
                     sizeof(double) * block_size * nvir * nvir * nE, ovvE_addr, &ovvE_addr);

        // Write out VOOE
        temp_ovoE->gemm(false, true, 1.0, block_iaQ, iEQ, 0.0);
        double** temp_ovoEp = temp_ovoE->pointer();
        for (size_t a = 0; a < nvir; a++){
            size_t local_i = 0;
            for (size_t i = bstart; i < bstart + block_size; i++){
                psio_address vooE_addr =
                    psio_get_address(PSIO_ZERO, sizeof(double) * (a * nocc + i) * nocc * nE);
                psio_->write(unit_, "EP2 I_vooE Integrals", (char*)temp_ovoEp[local_i * nvir + a],
                             sizeof(double) * nocc * nE, vooE_addr, &vooE_addr);
                local_i++;
            }
        }

    }
    // Blow away the temp tensors
    block_iaQ.reset();
    temp_ovvE.reset();
    temp_ovoE.reset();
    aEQ.reset();
    iEQ.reset();


    // ==> More Sizing <== /

    size_t aaE_size = memory_doubles_ / nvir * nvir * nE;
    if (aaE_size > nocc) aaE_size = nocc;
    // aaE_size = 2;

    size_t ooE_size = memory_doubles_ / nocc * nocc * nE;
    if (ooE_size > nvir) ooE_size = nvir;
    // ooE_size = 2;

    size_t aaE_nblocks = 1 + ((nocc - 1) / aaE_size);
    size_t ooE_nblocks = 1 + ((nvir - 1) / ooE_size);


    if (debug_ > 1){
        outfile->Printf("\n\n");
        outfile->Printf("aaE_size %zu\n", aaE_size);
        outfile->Printf("ooE_size    %zu\n", ooE_size);
        outfile->Printf("aaE_nblocks   %zu\n", aaE_nblocks);
        outfile->Printf("ooE_nblocks   %zu\n", ooE_nblocks);
        outfile->Printf("\n\n");
    }


    // ==> Iterate <== /
    outfile->Printf("  ==> Iterations <==\n\n");
    outfile->Printf( "   --------------------------------------------\n");
    outfile->Printf( "    Iter   Residual RMS      Max RMS    Remain\n");
    outfile->Printf( "   --------------------------------------------\n");

    std::vector<double> Eold(nE);
    std::vector<double> Esigma(nE);
    std::vector<double> Ederiv(nE);
    std::vector<double> Eerror(nE);

    // thread info
    std::vector<std::vector<double>> deriv_temps(num_threads_, std::vector<double>(nE));
    std::vector<std::vector<double>> sigma_temps(num_threads_, std::vector<double>(nE));
    size_t rank = 0;

    for (size_t iter = 0; iter < max_iter_; iter++) {

        // Reset data for loop
        for (size_t i = 0; i < nE; i++) {
            Esigma[i] = 0.0;
            Ederiv[i] = 0.0;
            Eold[i] = Enew[i];
            for (size_t j = 0; j < num_threads_; j++){
                deriv_temps[j][i] = 0.0;
                sigma_temps[j][i] = 0.0;
            }
        }


        // => Excitations <= //
        // sigma <= (Eabi - Ebai) * Eabi / (E - v - v + o)

        ovvE_addr = psio_get_address(PSIO_ZERO, 0);
        SharedMatrix I_ovvE(new Matrix(aaE_size * nvir, nvir * nE));
        double** I_ovvEp = I_ovvE->pointer();

        for (size_t i_block = 0; i_block < aaE_nblocks; i_block++){
            size_t i_start = aaE_size * i_block;
            size_t ib_size = aaE_size;
            if ((i_start + aaE_size) > nocc){
                ib_size = nocc - i_start;
            }

            psio_->read(unit_, "EP2 I_ovvE Integrals", (char*)I_ovvEp[0],
                        (sizeof(double) * ib_size * nvir * nvir * nE), ovvE_addr, &ovvE_addr);

            #pragma omp parallel for private(rank) schedule(dynamic,1) collapse(2) num_threads(num_threads_)
            for (size_t i = 0; i < ib_size; i++) {
                for (size_t a = 0; a < nvir; a++) {
                    #ifdef _OPENMP
                        rank = omp_get_thread_num();
                    #endif
                    # pragma omp simd collapse(2)
                    for (size_t b = 0; b < nvir; b++) {
                        for (size_t e = 0; e < nE; e++) {
                            double Eabi = I_ovvEp[i * nvir + b][a * nE + e];
                            double Ebai = I_ovvEp[i * nvir + a][b * nE + e];
                            double numer = (2.0 * Eabi - Ebai) * Eabi;
                            double denom = (denom_E[e] - eps_vir[a] - eps_vir[b] + eps_occ[i_start + i]);

                            sigma_temps[rank][e] += numer / denom;
                            deriv_temps[rank][e] += numer / (denom * denom);
                        }
                    }
                }
            }
        }
        I_ovvE.reset();

        // => De-excitations <= //
        // sigma <= (Eija - Ejia) * Eija / (E - o - o + v)

        vooE_addr = psio_get_address(PSIO_ZERO, 0);
        SharedMatrix I_vooE(new Matrix(ooE_size * nocc, nocc * nE));
        double** I_vooEp = I_vooE->pointer();

        for (size_t a_block = 0; a_block < ooE_nblocks; a_block++){
            size_t a_start = ooE_size * a_block;
            size_t ab_size = ooE_size;
            if ((a_start + ooE_size) > nvir){
                ab_size = nvir - a_start;
            }

            psio_->read(unit_, "EP2 I_vooE Integrals", (char*)I_vooEp[0],
                        sizeof(double) * ab_size * nocc * nocc * nE, vooE_addr, &vooE_addr);

            #pragma omp parallel for private(rank) schedule(dynamic,1) collapse(2) num_threads(num_threads_)
            for (size_t a = 0; a < ab_size; a++) {
                for (size_t i = 0; i < nocc; i++) {
                    #ifdef _OPENMP
                        rank = omp_get_thread_num();
                    #endif
                    # pragma omp simd collapse(2)
                    for (size_t j = 0; j < nocc; j++) {
                        for (size_t e = 0; e < nE; e++) {
                            double Eija = I_vooEp[a * nocc + j][i * nE + e];
                            double Ejia = I_vooEp[a * nocc + i][j * nE + e];
                            double numer = (2.0 * Eija - Ejia) * Eija;
                            double denom = (denom_E[e] - eps_occ[i] - eps_occ[j] + eps_vir[a_start + a]);

                            sigma_temps[rank][e] += numer / denom;
                            deriv_temps[rank][e] += numer / (denom * denom);
                        }
                    }
                }
            }
        }
        I_vooE.reset();

        // Sum up thread data
        for (size_t i = 0; i < nE; i++) {
            for (size_t j = 0; j < num_threads_; j++){
                Ederiv[i] += deriv_temps[j][i];
                Esigma[i] += sigma_temps[j][i];
            }
        }

        // Update
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

        }
        mean_error /= (double)nE;
        // printf("\n");

        outfile->Printf("    %3zu %14.8f %14.8f   %4zu\n", (iter + 1), mean_error, max_error, nremain);

        if (max_error < conv_thresh_) break;

        for (size_t i = 0; i < nE; i++) {
            // printf("%8.2f %8.2f |  ", Esigma[i], Ederiv[i]);
            // Reset E's
            Eold[i] = Enew[i];
            Esigma[i] = 0.0;
            Ederiv[i] = 0.0;
        }
        // printf("\n");
    }
    outfile->Printf( "   --------------------------------------------\n\n");
    psio->close(unit_, 0);

    // Build output array and remap symmetry
    std::vector<std::vector<std::pair<double, double>>> ret;
    for (size_t h = 0; h < nirrep_; h++) {
        std::vector<std::pair<double, double>> tmp(solve_orbs[h].size());
        ret.push_back(tmp);
    }

    for (size_t k = 0; k < orb_positions.size(); k++) {
        size_t h = std::get<1>(orb_positions[k]);
        size_t i = std::get<3>(orb_positions[k]);
        ret[h][i].first = Enew[k];
        ret[h][i].second = 1.0 / (1.0 - Ederiv[k]);
    }

    return ret;

}

}}
