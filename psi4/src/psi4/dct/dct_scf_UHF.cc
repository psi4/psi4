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

#include "dct.h"

#include <algorithm>
#include <cmath>
#include <map>

#include "psi4/libiwl/iwl.hpp"
#include "psi4/libdpd/dpd.h"
#include "psi4/libqt/qt.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/wavefunction.h"
#include "psi4/libmints/mintshelper.h"
#include "psi4/libtrans/integraltransform.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/liboptions/liboptions.h"

#include "psi4/psifiles.h"

namespace psi {
namespace dct {

/**
 * Checks to make sure that the phase, and ordering, of the MOs is consistent
 * with the previous set.
 *
 * @return Whether the phase correction was successful
 */
bool DCTSolver::correct_mo_phases(bool dieOnError) {
    dct_timer_on("DCTSolver::correct_mo_phases()");

    Matrix temp("temp", nirrep_, nsopi_, nmopi_);
    Matrix overlap("Old - New Overlap", nirrep_, nmopi_, nmopi_);

    bool error = correct_mo_phase_spincase(temp, overlap, *old_ca_, *Ca_, dieOnError);
    error = error && correct_mo_phase_spincase(temp, overlap, *old_cb_, *Cb_, dieOnError);

    dct_timer_off("DCTSolver::correct_mo_phases()");
    return error;
}

bool DCTSolver::correct_mo_phase_spincase(Matrix& temp, Matrix& overlap, const Matrix& old_C, Matrix& C, bool dieOnError) const {
    temp.gemm(false, false, 1.0, ao_s_, C, 0.0);
    overlap.gemm(true, false, 1.0, old_C, temp, 0.0);

    temp.copy(C);
    for (int h = 0; h < nirrep_; ++h) {
        std::map<int, int> mosUsed;
        for (int oldMO = 0; oldMO < nmopi_[h]; ++oldMO) {
            int bestMO = 0;
            double maximumProjection = 0.0;
            double prefactor = 0.0;
            for (int newMO = 0; newMO < nmopi_[h]; ++newMO) {
                double val = overlap.get(h, oldMO, newMO);
                if (std::fabs(val) > maximumProjection) {
                    maximumProjection = std::fabs(val);
                    bestMO = newMO;
                    prefactor = val < 0.0 ? -1.0 : 1.0;
                }
            }
            // Now we've found the MO to use, check it's not been used already then
            // copy it over.
            if (mosUsed[bestMO]++) {
                if (dieOnError) {
                    overlap.print();
                    old_C.print();
                    temp.print();
                    throw SanityCheckError("Two new MOs most resemble the same old MO. This has to be a bug.", __FILE__, __LINE__);
                } else {
                    // Copy Ca back from temp
                    C.copy(temp);
                    return false;
                }
            }
            for (int so = 0; so < nsopi_[h]; ++so) {
                C.set(h, so, oldMO, prefactor * temp.get(h, so, bestMO));
            }
        }
    }
    return true;
}

/**
 * Reads the orbitals and related quantities from the reference wavefunction
 * and reads the one-electron integrals from PSIO.
 */
void DCTSolver::initialize_orbitals_from_reference_U() {
    dct_timer_on("DCTSolver::scf_guess");

    epsilon_a_->copy(*reference_wavefunction_->epsilon_a());
    epsilon_b_->copy(*reference_wavefunction_->epsilon_b());
    Ca_->copy(reference_wavefunction_->Ca());
    Cb_->copy(reference_wavefunction_->Cb());
    moFa_->copy(reference_wavefunction_->Fa());
    moFa_->transform(Ca_);
    moFb_->copy(reference_wavefunction_->Fb());
    moFb_->transform(Cb_);
    update_scf_density();

    dct_timer_off("DCTSolver::scf_guess");
}

/**
 * Computes the SCF energy from the latest Fock and density matrices.
 * WARNING! This quantity is a misnomer from earlier days of the theory.
 * "SCF" here means "excluding 2RDM cumulant, including 1RDM and 1RDM products."
 */
void DCTSolver::compute_scf_energy() {
    dct_timer_on("DCTSolver::compute_scf_energy");

    // Escf = eNuc + 0.5 * (H + F) * (kappa + tau)
    scf_energy_ = enuc_;
    scf_energy_ += 0.5 * kappa_so_a_->vector_dot(so_h_);
    scf_energy_ += 0.5 * kappa_so_b_->vector_dot(so_h_);
    scf_energy_ += 0.5 * tau_so_a_->vector_dot(so_h_);
    scf_energy_ += 0.5 * tau_so_b_->vector_dot(so_h_);

    if (options_.get_str("DCT_TYPE") == "DF" && options_.get_str("AO_BASIS") == "NONE") {
        scf_energy_ += 0.5 * mo_gammaA_.vector_dot(moFa_);
        scf_energy_ += 0.5 * mo_gammaB_.vector_dot(moFb_);
    } else {
        scf_energy_ += 0.5 * kappa_so_a_->vector_dot(Fa_);
        scf_energy_ += 0.5 * kappa_so_b_->vector_dot(Fb_);
        scf_energy_ += 0.5 * tau_so_a_->vector_dot(Fa_);
        scf_energy_ += 0.5 * tau_so_b_->vector_dot(Fb_);
    }

    dct_timer_off("DCTSolver::compute_scf_energy");
}

/**
 * Computes the SCF error vector by transforming the Fock matrices to the
 * MO basis and computing [F, Kappa], and the RMS value of this quantity.
 * @return RMS error
 */
double DCTSolver::compute_scf_error_vector() {
    dct_timer_on("DCTSolver::compute_scf_error_vector");

    size_t nElements = 0;
    double sumOfSquares = 0.0;
    auto tmp1 = Matrix("tmp1", nirrep_, nsopi_, nsopi_);
    auto tmp2 = Matrix("tmp2", nirrep_, nsopi_, nsopi_);
    // form FDS
    tmp1.gemm(false, false, 1.0, kappa_so_a_, ao_s_, 0.0);
    scf_error_a_->gemm(false, false, 1.0, Fa_, tmp1, 0.0);
    // form SDF
    tmp1.gemm(false, false, 1.0, kappa_so_a_, Fa_, 0.0);
    tmp2.gemm(false, false, 1.0, ao_s_, tmp1, 0.0);
    scf_error_a_->subtract(tmp2);
    // Orthogonalize
    scf_error_a_->transform(s_half_inv_);

    // form FDS
    tmp1.gemm(false, false, 1.0, kappa_so_b_, ao_s_, 0.0);
    scf_error_b_->gemm(false, false, 1.0, Fb_, tmp1, 0.0);
    // form SDF
    tmp1.gemm(false, false, 1.0, kappa_so_b_, Fb_, 0.0);
    tmp2.gemm(false, false, 1.0, ao_s_, tmp1, 0.0);
    scf_error_b_->subtract(tmp2);
    // Orthogonalize
    scf_error_b_->transform(s_half_inv_);

    for (int h = 0; h < nirrep_; ++h) {
        for (int p = 0; p < nsopi_[h]; ++p) {
            for (int q = 0; q < nsopi_[h]; ++q) {
                nElements += 2;
                sumOfSquares += pow(scf_error_a_->get(h, p, q), 2.0);
                sumOfSquares += pow(scf_error_b_->get(h, p, q), 2.0);
            }
        }
    }
    dct_timer_off("DCTSolver::compute_scf_error_vector");
    return sqrt(sumOfSquares / nElements);
}

/**
 * Uses the MO coefficients to form the SCF density matrices in the SO basis.
 * @param Whether to damp the update or not
 * @return RMS density change
 */
double DCTSolver::update_scf_density(bool damp) {
    dct_timer_on("DCTSolver::update_scf_density");

    double dampingFactor = options_.get_double("DAMPING_PERCENTAGE");
    double newFraction = damp ? 1.0 : 1.0 - dampingFactor / 100.0;
    size_t nElements = 0;
    double sumOfSquares = 0.0;
    Matrix old(kappa_so_a_);
    for (int h = 0; h < nirrep_; ++h) {
        for (int mu = 0; mu < nsopi_[h]; ++mu) {
            for (int nu = 0; nu < nsopi_[h]; ++nu) {
                double val = 0.0;
                for (int i = 0; i < naoccpi_[h]; ++i) val += Ca_->get(h, mu, i) * Ca_->get(h, nu, i);
                kappa_so_a_->set(h, mu, nu, newFraction * val + (1.0 - newFraction) * kappa_so_a_->get(h, mu, nu));
                ++nElements;
                sumOfSquares += pow(val - old.get(h, mu, nu), 2.0);
            }
        }
    }
    old.copy(kappa_so_b_);
    for (int h = 0; h < nirrep_; ++h) {
        for (int mu = 0; mu < nsopi_[h]; ++mu) {
            for (int nu = 0; nu < nsopi_[h]; ++nu) {
                double val = 0.0;
                for (int i = 0; i < nboccpi_[h]; ++i) val += Cb_->get(h, mu, i) * Cb_->get(h, nu, i);
                kappa_so_b_->set(h, mu, nu, newFraction * val + (1.0 - newFraction) * kappa_so_b_->get(h, mu, nu));
                ++nElements;
                sumOfSquares += pow(val - old.get(h, mu, nu), 2.0);
            }
        }
    }
    // We're not converged until the RMS error vector *and* the RMS density
    // changes are below the threshold
    dct_timer_off("DCTSolver::update_scf_density");

    return sqrt(sumOfSquares / nElements);
}

/**
 * Builds the G matrix for Fock matrix, the external potential (tau contracted with the integrals)
 * and other tensors, if requested, out-of-core using the SO integrals.
 * Also builds the AO basis tensors for the <VV||VV>, <vv||vv>, and <Vv|Vv> terms in the G and X intermediates
 * All quantities are built simultaneously to reduce I/O.
 */
void DCTSolver::process_so_ints() {
    dct_timer_on("DCTSolver::process_so_ints");

    IWL *iwl = new IWL(psio_.get(), PSIF_SO_TEI, int_tolerance_, 1, 1);

    Label *lblptr = iwl->labels();
    Value *valptr = iwl->values();

    double *Da = init_array(ntriso_);
    double *Db = init_array(ntriso_);
    double *Ga = init_array(ntriso_);
    double *Gb = init_array(ntriso_);

    auto opdm_a = kappa_so_a_->clone();
    opdm_a->add(tau_so_a_);
    auto opdm_b = kappa_so_b_->clone();
    opdm_b->add(tau_so_b_);
    int soOffset = 0;
    for (int h = 0; h < nirrep_; ++h) {
        for (int mu = 0; mu < nsopi_[h]; ++mu) {
            for (int nu = 0; nu <= mu; ++nu) {
                int muNu = INDEX((nu + soOffset), (mu + soOffset));
                Da[muNu] = opdm_a->get(h, mu, nu);
                Db[muNu] = opdm_b->get(h, mu, nu);
            }
        }
        soOffset += nsopi_[h];
    }

    double value;
    int Gc, Gd;
    int pqArr, qpArr, rsArr, srArr, qrArr, rqArr;
    int qsArr, sqArr, psArr, spArr, prArr, rpArr;
    int offset, labelIndex, p, q, r, s, h, counter;
    int **pq_row_start, **CD_row_start, **Cd_row_start, **cd_row_start;
    dpdbuf4 tau_temp, lambda;
    dpdbuf4 tau1_AO_aa, tau2_AO_aa;
    dpdbuf4 tau1_AO_ab, tau2_AO_ab;
    dpdbuf4 tau1_AO_bb, tau2_AO_bb;

    bool buildTensors = (options_.get_str("AO_BASIS") == "DISK");

    if (buildTensors) {
        counter = 0;

        // Build the offset arrays needed for the DGEMM in half_transform
        pq_row_start = init_int_matrix(nirrep_, nirrep_);
        CD_row_start = init_int_matrix(nirrep_, nirrep_);
        cd_row_start = init_int_matrix(nirrep_, nirrep_);
        Cd_row_start = init_int_matrix(nirrep_, nirrep_);
        for (h = 0; h < nirrep_; ++h) {
            for (Gc = 0, offset = 0; Gc < nirrep_; ++Gc) {
                Gd = Gc ^ h;
                pq_row_start[h][Gc] = offset;
                offset += nsopi_[Gc] * nsopi_[Gd];
            }
            for (Gc = 0, offset = 0; Gc < nirrep_; ++Gc) {
                Gd = Gc ^ h;
                CD_row_start[h][Gc] = offset;
                offset += navirpi_[Gc] * navirpi_[Gd];
            }
            for (Gc = 0, offset = 0; Gc < nirrep_; ++Gc) {
                Gd = Gc ^ h;
                Cd_row_start[h][Gc] = offset;
                offset += navirpi_[Gc] * nbvirpi_[Gd];
            }
            for (Gc = 0, offset = 0; Gc < nirrep_; ++Gc) {
                Gd = Gc ^ h;
                cd_row_start[h][Gc] = offset;
                offset += nbvirpi_[Gc] * nbvirpi_[Gd];
            }
        }

        dpd_set_default(_ints->get_dpd_id());

        /********** AA ***********/
        global_dpd_->buf4_init(&lambda, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O>O]-"), ID("[V>V]-"), 0,
                               "Amplitude <OO|VV>");
        global_dpd_->buf4_init(&tau1_AO_aa, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[n,n]"), ID("[O,O]"), ID("[n,n]"), 0,
                               "tau1AO <OO|nn>");
        global_dpd_->buf4_scm(&tau1_AO_aa, 0.0);
        half_transform(&tau1_AO_aa, &lambda, avir_c_, avir_c_, navirpi_, navirpi_, pq_row_start, CD_row_start, true,
                       1.0, 0.0);
        global_dpd_->buf4_close(&lambda);
        global_dpd_->buf4_close(&tau1_AO_aa);

        // Now sort for better memory access patterns
        global_dpd_->buf4_init(&tau1_AO_aa, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[n,n]"), ID("[O,O]"), ID("[n,n]"), 1,
                               "tau1AO <OO|nn>");
        global_dpd_->buf4_sort(&tau1_AO_aa, PSIF_DCT_DPD, rspq, ID("[n,n]"), ID("[O,O]"), "tau1AO <nn|OO>");

        global_dpd_->buf4_close(&tau1_AO_aa);

        // Now reopen the two AO dpd_buf4's
        global_dpd_->buf4_init(&tau1_AO_aa, PSIF_DCT_DPD, 0, ID("[n,n]"), ID("[O,O]"), ID("[n,n]"), ID("[O,O]"), 0,
                               "tau1AO <nn|OO>");
        global_dpd_->buf4_init(&tau2_AO_aa, PSIF_DCT_DPD, 0, ID("[n,n]"), ID("[O,O]"), ID("[n,n]"), ID("[O,O]"), 0,
                               "tau2AO <nn|OO>");
        global_dpd_->buf4_scm(&tau2_AO_aa, 0.0);

        /********** BB ***********/
        global_dpd_->buf4_init(&lambda, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o>o]-"), ID("[v>v]-"), 0,
                               "Amplitude <oo|vv>");
        global_dpd_->buf4_init(&tau1_AO_bb, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[n,n]"), ID("[o,o]"), ID("[n,n]"), 0,
                               "tau1AO <oo|nn>");
        global_dpd_->buf4_scm(&tau1_AO_bb, 0.0);
        half_transform(&tau1_AO_bb, &lambda, bvir_c_, bvir_c_, nbvirpi_, nbvirpi_, pq_row_start, cd_row_start, true,
                       1.0, 0.0);
        global_dpd_->buf4_close(&lambda);
        global_dpd_->buf4_close(&tau1_AO_bb);

        // Now sort for better memory access patterns
        global_dpd_->buf4_init(&tau1_AO_bb, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[n,n]"), ID("[o,o]"), ID("[n,n]"), 1,
                               "tau1AO <oo|nn>");
        global_dpd_->buf4_sort(&tau1_AO_bb, PSIF_DCT_DPD, rspq, ID("[n,n]"), ID("[o,o]"), "tau1AO <nn|oo>");
        global_dpd_->buf4_close(&tau1_AO_bb);

        // Now reopen the two AO dpd_buf4's
        global_dpd_->buf4_init(&tau1_AO_bb, PSIF_DCT_DPD, 0, ID("[n,n]"), ID("[o,o]"), ID("[n,n]"), ID("[o,o]"), 0,
                               "tau1AO <nn|oo>");
        global_dpd_->buf4_init(&tau2_AO_bb, PSIF_DCT_DPD, 0, ID("[n,n]"), ID("[o,o]"), ID("[n,n]"), ID("[o,o]"), 0,
                               "tau2AO <nn|oo>");
        global_dpd_->buf4_scm(&tau2_AO_bb, 0.0);

        /********** AB ***********/
        global_dpd_->buf4_init(&lambda, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0,
                               "Amplitude <Oo|Vv>");
        global_dpd_->buf4_init(&tau1_AO_ab, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[n,n]"), ID("[O,o]"), ID("[n,n]"), 0,
                               "tau1AO <Oo|nn>");
        global_dpd_->buf4_scm(&tau1_AO_ab, 0.0);
        half_transform(&tau1_AO_ab, &lambda, avir_c_, bvir_c_, navirpi_, nbvirpi_, pq_row_start, Cd_row_start, true,
                       1.0, 0.0);
        global_dpd_->buf4_close(&lambda);
        global_dpd_->buf4_close(&tau1_AO_ab);

        // Now sort for better memory access patterns
        global_dpd_->buf4_init(&tau1_AO_ab, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[n,n]"), ID("[O,o]"), ID("[n,n]"), 0,
                               "tau1AO <Oo|nn>");
        global_dpd_->buf4_sort(&tau1_AO_ab, PSIF_DCT_DPD, rspq, ID("[n,n]"), ID("[O,o]"), "tau1AO <nn|Oo>");
        global_dpd_->buf4_close(&tau1_AO_ab);

        // Reopen the two AO dpd_buf4's
        global_dpd_->buf4_init(&tau1_AO_ab, PSIF_DCT_DPD, 0, ID("[n,n]"), ID("[O,o]"), ID("[n,n]"), ID("[O,o]"), 0,
                               "tau1AO <nn|Oo>");
        global_dpd_->buf4_init(&tau2_AO_ab, PSIF_DCT_DPD, 0, ID("[n,n]"), ID("[O,o]"), ID("[n,n]"), ID("[O,o]"), 0,
                               "tau2AO <nn|Oo>");
        global_dpd_->buf4_scm(&tau2_AO_ab, 0.0);

        // Now put stuff in memory
        for (int h = 0; h < nirrep_; ++h) {
            global_dpd_->buf4_mat_irrep_init(&tau1_AO_aa, h);
            global_dpd_->buf4_mat_irrep_rd(&tau1_AO_aa, h);
            global_dpd_->buf4_mat_irrep_init(&tau2_AO_aa, h);

            global_dpd_->buf4_mat_irrep_init(&tau1_AO_bb, h);
            global_dpd_->buf4_mat_irrep_rd(&tau1_AO_bb, h);
            global_dpd_->buf4_mat_irrep_init(&tau2_AO_bb, h);

            global_dpd_->buf4_mat_irrep_init(&tau1_AO_ab, h);
            global_dpd_->buf4_mat_irrep_rd(&tau1_AO_ab, h);
            global_dpd_->buf4_mat_irrep_init(&tau2_AO_ab, h);
        }
    }

    bool lastBuffer;
    do {
        lastBuffer = iwl->last_buffer();
        for (int index = 0; index < iwl->buffer_count(); ++index) {
            labelIndex = 4 * index;
            p = std::abs((int)lblptr[labelIndex++]);
            q = (int)lblptr[labelIndex++];
            r = (int)lblptr[labelIndex++];
            s = (int)lblptr[labelIndex++];
            value = (double)valptr[index];
            if (buildTensors) {
                AO_contribute(&tau1_AO_aa, &tau2_AO_aa, p, q, r, s, value);
                AO_contribute(&tau1_AO_bb, &tau2_AO_bb, p, q, r, s, value);
                AO_contribute(&tau1_AO_ab, &tau2_AO_ab, p, q, r, s, value);
                ++counter;
            }

            qpArr = pqArr = INDEX(p, q);
            srArr = rsArr = INDEX(r, s);
            prArr = rpArr = INDEX(p, r);
            qsArr = sqArr = INDEX(q, s);
            spArr = psArr = INDEX(p, s);
            qrArr = rqArr = INDEX(q, r);

            /* (pq|rs) */
            Ga[rsArr] += (Da[pqArr] + Db[pqArr]) * value;
            Gb[rsArr] += (Da[pqArr] + Db[pqArr]) * value;
            if (q >= r) {
                Ga[qrArr] -= Da[psArr] * value;
                Gb[qrArr] -= Db[psArr] * value;
            }

            if (p != q && r != s && pqArr != rsArr) {
                /* (pq|sr) */
                if (s >= r) {
                    Ga[srArr] += (Da[pqArr] + Db[pqArr]) * value;
                    Gb[srArr] += (Da[pqArr] + Db[pqArr]) * value;
                }
                if (q >= s) {
                    Ga[qsArr] -= Da[prArr] * value;
                    Gb[qsArr] -= Db[prArr] * value;
                }

                /* (qp|rs) */
                if (r >= s) {
                    Ga[rsArr] += (Da[qpArr] + Db[qpArr]) * value;
                    Gb[rsArr] += (Da[qpArr] + Db[qpArr]) * value;
                }
                if (p >= r) {
                    Ga[prArr] -= Da[qsArr] * value;
                    Gb[prArr] -= Db[qsArr] * value;
                }

                /* (qp|sr) */
                if (s >= r) {
                    Ga[srArr] += (Da[qpArr] + Db[qpArr]) * value;
                    Gb[srArr] += (Da[qpArr] + Db[qpArr]) * value;
                }
                if (p >= s) {
                    Ga[psArr] -= Da[qrArr] * value;
                    Gb[psArr] -= Db[qrArr] * value;
                }

                /* (rs|pq) */
                if (p >= q) {
                    Ga[pqArr] += (Da[rsArr] + Db[rsArr]) * value;
                    Gb[pqArr] += (Da[rsArr] + Db[rsArr]) * value;
                }
                if (s >= p) {
                    Ga[spArr] -= Da[rqArr] * value;
                    Gb[spArr] -= Db[rqArr] * value;
                }

                /* (sr|pq) */
                if (p >= q) {
                    Ga[pqArr] += (Da[srArr] + Db[srArr]) * value;
                    Gb[pqArr] += (Da[srArr] + Db[srArr]) * value;
                }
                if (r >= p) {
                    Ga[rpArr] -= Da[sqArr] * value;
                    Gb[rpArr] -= Db[sqArr] * value;
                }

                /* (rs|qp) */
                if (q >= p) {
                    Ga[qpArr] += (Da[rsArr] + Db[rsArr]) * value;
                    Gb[qpArr] += (Da[rsArr] + Db[rsArr]) * value;
                }
                if (s >= q) {
                    Ga[sqArr] -= Da[rpArr] * value;
                    Gb[sqArr] -= Db[rpArr] * value;
                }

                /* (sr|qp) */
                if (q >= p) {
                    Ga[qpArr] += (Da[srArr] + Db[srArr]) * value;
                    Gb[qpArr] += (Da[srArr] + Db[srArr]) * value;
                }
                if (r >= q) {
                    Ga[rqArr] -= Da[spArr] * value;
                    Gb[rqArr] -= Db[spArr] * value;
                }
            } else if (p != q && r != s && pqArr == rsArr) {
                /* (pq|sr) */
                if (s >= r) {
                    Ga[srArr] += (Da[pqArr] + Db[pqArr]) * value;
                    Gb[srArr] += (Da[pqArr] + Db[pqArr]) * value;
                }
                if (q >= s) {
                    Ga[qsArr] -= Da[prArr] * value;
                    Gb[qsArr] -= Db[prArr] * value;
                }
                /* (qp|rs) */
                if (r >= s) {
                    Ga[rsArr] += (Da[qpArr] + Db[qpArr]) * value;
                    Gb[rsArr] += (Da[qpArr] + Db[qpArr]) * value;
                }
                if (p >= r) {
                    Ga[prArr] -= Da[qsArr] * value;
                    Gb[prArr] -= Db[qsArr] * value;
                }

                /* (qp|sr) */
                if (s >= r) {
                    Ga[srArr] += (Da[qpArr] + Db[qpArr]) * value;
                    Gb[srArr] += (Da[qpArr] + Db[qpArr]) * value;
                }
                if (p >= s) {
                    Ga[psArr] -= Da[qrArr] * value;
                    Gb[psArr] -= Db[qrArr] * value;
                }
            } else if (p != q && r == s) {
                /* (qp|rs) */
                if (r >= s) {
                    Ga[rsArr] += (Da[qpArr] + Db[qpArr]) * value;
                    Gb[rsArr] += (Da[qpArr] + Db[qpArr]) * value;
                }
                if (p >= r) {
                    Ga[prArr] -= Da[qsArr] * value;
                    Gb[prArr] -= Db[qsArr] * value;
                }

                /* (rs|pq) */
                if (p >= q) {
                    Ga[pqArr] += (Da[rsArr] + Db[rsArr]) * value;
                    Gb[pqArr] += (Da[rsArr] + Db[rsArr]) * value;
                }
                if (s >= p) {
                    Ga[spArr] -= Da[rqArr] * value;
                    Gb[spArr] -= Db[rqArr] * value;
                }

                /* (rs|qp) */
                if (q >= p) {
                    Ga[qpArr] += (Da[rsArr] + Db[rsArr]) * value;
                    Gb[qpArr] += (Da[rsArr] + Db[rsArr]) * value;
                }
                if (s >= q) {
                    Ga[sqArr] -= Da[rpArr] * value;
                    Gb[sqArr] -= Db[rpArr] * value;
                }
            } else if (p == q && r != s) {
                /* (pq|sr) */
                if (s >= r) {
                    Ga[srArr] += (Da[pqArr] + Db[pqArr]) * value;
                    Gb[srArr] += (Da[pqArr] + Db[pqArr]) * value;
                }
                if (q >= s) {
                    Ga[qsArr] -= Da[prArr] * value;
                    Gb[qsArr] -= Db[prArr] * value;
                }

                /* (rs|pq) */
                if (p >= q) {
                    Ga[pqArr] += (Da[rsArr] + Db[rsArr]) * value;
                    Gb[pqArr] += (Da[rsArr] + Db[rsArr]) * value;
                }
                if (s >= p) {
                    Ga[spArr] -= Da[rqArr] * value;
                    Gb[spArr] -= Db[rqArr] * value;
                }

                /* (sr|pq) */
                if (p >= q) {
                    Ga[pqArr] += (Da[srArr] + Db[srArr]) * value;
                    Gb[pqArr] += (Da[srArr] + Db[srArr]) * value;
                }
                if (r >= p) {
                    Ga[rpArr] -= Da[sqArr] * value;
                    Gb[rpArr] -= Db[sqArr] * value;
                }
            } else if (p == q && r == s && pqArr != rsArr) {
                /* (rs|pq) */
                if (p >= q) {
                    Ga[pqArr] += (Da[rsArr] + Db[rsArr]) * value;
                    Gb[pqArr] += (Da[rsArr] + Db[rsArr]) * value;
                }
                if (s >= p) {
                    Ga[spArr] -= Da[rqArr] * value;
                    Gb[spArr] -= Db[rqArr] * value;
                }
            }
        } /* end loop through current buffer */
        if (!lastBuffer) iwl->fetch();
    } while (!lastBuffer);
    iwl->set_keep_flag(true);
    delete iwl;
    if (buildTensors) {
        if (print_ > 1) {
            outfile->Printf("Processed %d SO integrals each for AA, BB, and AB\n", counter);
        }
        for (int h = 0; h < nirrep_; ++h) {
            global_dpd_->buf4_mat_irrep_wrt(&tau2_AO_aa, h);
            global_dpd_->buf4_mat_irrep_close(&tau1_AO_aa, h);
            global_dpd_->buf4_mat_irrep_close(&tau2_AO_aa, h);

            global_dpd_->buf4_mat_irrep_wrt(&tau2_AO_bb, h);
            global_dpd_->buf4_mat_irrep_close(&tau1_AO_bb, h);
            global_dpd_->buf4_mat_irrep_close(&tau2_AO_bb, h);

            global_dpd_->buf4_mat_irrep_wrt(&tau2_AO_ab, h);
            global_dpd_->buf4_mat_irrep_close(&tau1_AO_ab, h);
            global_dpd_->buf4_mat_irrep_close(&tau2_AO_ab, h);
        }

        global_dpd_->buf4_close(&tau1_AO_aa);
        global_dpd_->buf4_close(&tau1_AO_bb);
        global_dpd_->buf4_close(&tau1_AO_ab);
        global_dpd_->buf4_close(&tau2_AO_aa);
        global_dpd_->buf4_close(&tau2_AO_bb);
        global_dpd_->buf4_close(&tau2_AO_ab);

        /********** AA ***********/
        global_dpd_->buf4_init(&tau2_AO_aa, PSIF_DCT_DPD, 0, ID("[n,n]"), ID("[O,O]"), ID("[n,n]"), ID("[O,O]"), 0,
                               "tau2AO <nn|OO>");
        global_dpd_->buf4_sort(&tau2_AO_aa, PSIF_DCT_DPD, rspq, ID("[O,O]"), ID("[n,n]"), "tau2AO <OO|nn>");
        global_dpd_->buf4_close(&tau2_AO_aa);
        global_dpd_->buf4_init(&tau2_AO_aa, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[n,n]"), ID("[O,O]"), ID("[n,n]"), 0,
                               "tau2AO <OO|nn>");
        global_dpd_->buf4_init(&tau_temp, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                               "tau(temp) <OO|VV>");
        global_dpd_->buf4_scm(&tau_temp, 0.0);
        half_transform(&tau2_AO_aa, &tau_temp, avir_c_, avir_c_, navirpi_, navirpi_, pq_row_start, CD_row_start, false,
                       0.5, 0.0);
        global_dpd_->buf4_close(&tau2_AO_aa);
        global_dpd_->buf4_close(&tau_temp);

        /********** BB ***********/
        global_dpd_->buf4_init(&tau2_AO_bb, PSIF_DCT_DPD, 0, ID("[n,n]"), ID("[o,o]"), ID("[n,n]"), ID("[o,o]"), 0,
                               "tau2AO <nn|oo>");
        global_dpd_->buf4_sort(&tau2_AO_bb, PSIF_DCT_DPD, rspq, ID("[o,o]"), ID("[n,n]"), "tau2AO <oo|nn>");
        global_dpd_->buf4_close(&tau2_AO_bb);
        global_dpd_->buf4_init(&tau2_AO_bb, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[n,n]"), ID("[o,o]"), ID("[n,n]"), 0,
                               "tau2AO <oo|nn>");
        global_dpd_->buf4_init(&tau_temp, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o,o]"), ID("[v,v]"), 0,
                               "tau(temp) <oo|vv>");
        global_dpd_->buf4_scm(&tau_temp, 0.0);
        half_transform(&tau2_AO_bb, &tau_temp, bvir_c_, bvir_c_, nbvirpi_, nbvirpi_, pq_row_start, cd_row_start, false,
                       0.5, 0.0);
        global_dpd_->buf4_close(&tau2_AO_bb);
        global_dpd_->buf4_close(&tau_temp);

        /********** AB ***********/
        global_dpd_->buf4_init(&tau2_AO_ab, PSIF_DCT_DPD, 0, ID("[n,n]"), ID("[O,o]"), ID("[n,n]"), ID("[O,o]"), 0,
                               "tau2AO <nn|Oo>");
        global_dpd_->buf4_sort(&tau2_AO_ab, PSIF_DCT_DPD, rspq, ID("[O,o]"), ID("[n,n]"), "tau2AO <Oo|nn>");
        global_dpd_->buf4_close(&tau2_AO_ab);
        global_dpd_->buf4_init(&tau2_AO_ab, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[n,n]"), ID("[O,o]"), ID("[n,n]"), 0,
                               "tau2AO <Oo|nn>");
        global_dpd_->buf4_init(&tau_temp, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0,
                               "tau(temp) <Oo|Vv>");
        global_dpd_->buf4_scm(&tau_temp, 0.0);
        half_transform(&tau2_AO_ab, &tau_temp, avir_c_, bvir_c_, navirpi_, nbvirpi_, pq_row_start, Cd_row_start, false,
                       1.0, 0.0);
        global_dpd_->buf4_close(&tau2_AO_ab);
        global_dpd_->buf4_close(&tau_temp);

        free_int_matrix(pq_row_start);
        free_int_matrix(CD_row_start);
        free_int_matrix(cd_row_start);
        free_int_matrix(Cd_row_start);
    }

    // Build the Fock matrices from the H and G matrices
    soOffset = 0;
    for (int h = 0; h < nirrep_; ++h) {
        for (int mu = 0; mu < nsopi_[h]; ++mu) {
            for (int nu = 0; nu <= mu; ++nu) {
                int muNu = INDEX((nu + soOffset), (mu + soOffset));
                double aVal = Ga[muNu];
                double bVal = Gb[muNu];
                Fa_->add(h, mu, nu, aVal);
                Fb_->add(h, mu, nu, bVal);
                if (mu != nu) {
                    Fa_->add(h, nu, mu, aVal);
                    Fb_->add(h, nu, mu, bVal);
                }
            }
        }
        soOffset += nsopi_[h];
    }

    free(Da);
    free(Db);
    free(Ga);
    free(Gb);

    dct_timer_off("DCTSolver::process_so_ints");
}

/*
 * Update the Fock operator, defined as h^p_q + gbar^pr_qs 1RDM^s_r
 * This is an important intermediate in -06 DCT theories.
 * We also construct the matrix's associated denominators. These appear in first-order update steps.
 */
void DCTSolver::update_fock() {
    dct_timer_on("DCTSolver::update_fock");

    dpdfile2 Gtau;

    moFa_->copy(so_h_);
    moFa_->transform(Ca_);

    moFb_->copy(so_h_);
    moFb_->transform(Cb_);

    // We already have the gbar * tau contraction computed, so let's just add it in.
    auto moG_tau = Matrix("GGamma in the MO basis", nirrep_, nmopi_, nmopi_);

    // Alpha occupied
    global_dpd_->file2_init(&Gtau, PSIF_DCT_DPD, 0, ID('O'), ID('O'), "GGamma <O|O>");
    moG_tau.set_block(slices_.at("ACTIVE_OCC_A"), Matrix(&Gtau));
    global_dpd_->file2_close(&Gtau);

    // Alpha virtual
    global_dpd_->file2_init(&Gtau, PSIF_DCT_DPD, 0, ID('V'), ID('V'), "GGamma <V|V>");
    moG_tau.set_block(slices_.at("ACTIVE_VIR_A"), Matrix(&Gtau));
    global_dpd_->file2_close(&Gtau);

    moFa_->add(moG_tau);

    // Beta occupied
    moG_tau.zero();
    global_dpd_->file2_init(&Gtau, PSIF_DCT_DPD, 0, ID('o'), ID('o'), "GGamma <o|o>");
    moG_tau.set_block(slices_.at("ACTIVE_OCC_B"), Matrix(&Gtau));
    global_dpd_->file2_close(&Gtau);

    // Beta virtual
    global_dpd_->file2_init(&Gtau, PSIF_DCT_DPD, 0, ID('v'), ID('v'), "GGamma <v|v>");
    moG_tau.set_block(slices_.at("ACTIVE_VIR_B"), Matrix(&Gtau));
    global_dpd_->file2_close(&Gtau);

    moFb_->add(moG_tau);

    // Write the MO basis Fock operator to the DPD file and update the denominators
    build_denominators();

    dct_timer_off("DCTSolver::update_fock");
}

/**
 * Builds the AO basis tensors for the <VV||VV>, <vv||vv>, and <Vv|Vv> terms in the G and X intermediates
 */
void DCTSolver::build_AO_tensors() {
    dct_timer_on("DCTSolver::build_tensors");

    IWL *iwl = new IWL(psio_.get(), PSIF_SO_TEI, int_tolerance_, 1, 1);

    Label *lblptr = iwl->labels();
    Value *valptr = iwl->values();

    double value;
    int Gc, Gd;
    int offset, labelIndex, p, q, r, s, h, counter;
    int **pq_row_start, **CD_row_start, **Cd_row_start, **cd_row_start;
    dpdbuf4 tau_temp, lambda;
    dpdbuf4 tau1_AO_aa, tau2_AO_aa;
    dpdbuf4 tau1_AO_ab, tau2_AO_ab;
    dpdbuf4 tau1_AO_bb, tau2_AO_bb;
    dpdfile2 s_aa_1, s_aa_2, s_aa_3, s_aa_4, tau;
    dpdfile2 s_bb_1, s_bb_2, s_bb_3, s_bb_4;

    counter = 0;

    // Build the offset arrays needed for the DGEMM in half_transform
    pq_row_start = init_int_matrix(nirrep_, nirrep_);
    CD_row_start = init_int_matrix(nirrep_, nirrep_);
    cd_row_start = init_int_matrix(nirrep_, nirrep_);
    Cd_row_start = init_int_matrix(nirrep_, nirrep_);
    for (h = 0; h < nirrep_; ++h) {
        for (Gc = 0, offset = 0; Gc < nirrep_; ++Gc) {
            Gd = Gc ^ h;
            pq_row_start[h][Gc] = offset;
            offset += nsopi_[Gc] * nsopi_[Gd];
        }
        for (Gc = 0, offset = 0; Gc < nirrep_; ++Gc) {
            Gd = Gc ^ h;
            CD_row_start[h][Gc] = offset;
            offset += navirpi_[Gc] * navirpi_[Gd];
        }
        for (Gc = 0, offset = 0; Gc < nirrep_; ++Gc) {
            Gd = Gc ^ h;
            Cd_row_start[h][Gc] = offset;
            offset += navirpi_[Gc] * nbvirpi_[Gd];
        }
        for (Gc = 0, offset = 0; Gc < nirrep_; ++Gc) {
            Gd = Gc ^ h;
            cd_row_start[h][Gc] = offset;
            offset += nbvirpi_[Gc] * nbvirpi_[Gd];
        }
    }

    dpd_set_default(_ints->get_dpd_id());

    /********** AA ***********/
    global_dpd_->buf4_init(&lambda, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O>O]-"), ID("[V>V]-"), 0,
                           "Amplitude <OO|VV>");
    global_dpd_->buf4_init(&tau1_AO_aa, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[n,n]"), ID("[O,O]"), ID("[n,n]"), 0,
                           "tau1AO <OO|nn>");
    global_dpd_->buf4_scm(&tau1_AO_aa, 0.0);
    half_transform(&tau1_AO_aa, &lambda, avir_c_, avir_c_, navirpi_, navirpi_, pq_row_start, CD_row_start, true, 1.0,
                   0.0);
    global_dpd_->buf4_close(&lambda);
    global_dpd_->buf4_close(&tau1_AO_aa);

    // Now sort for better memory access patterns
    global_dpd_->buf4_init(&tau1_AO_aa, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[n,n]"), ID("[O,O]"), ID("[n,n]"), 1,
                           "tau1AO <OO|nn>");
    global_dpd_->buf4_sort(&tau1_AO_aa, PSIF_DCT_DPD, rspq, ID("[n,n]"), ID("[O,O]"), "tau1AO <nn|OO>");

    global_dpd_->buf4_close(&tau1_AO_aa);

    // Now reopen the two AO dpd_buf4's
    global_dpd_->buf4_init(&tau1_AO_aa, PSIF_DCT_DPD, 0, ID("[n,n]"), ID("[O,O]"), ID("[n,n]"), ID("[O,O]"), 0,
                           "tau1AO <nn|OO>");
    global_dpd_->buf4_init(&tau2_AO_aa, PSIF_DCT_DPD, 0, ID("[n,n]"), ID("[O,O]"), ID("[n,n]"), ID("[O,O]"), 0,
                           "tau2AO <nn|OO>");
    global_dpd_->buf4_scm(&tau2_AO_aa, 0.0);

    /********** BB ***********/
    global_dpd_->buf4_init(&lambda, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o>o]-"), ID("[v>v]-"), 0,
                           "Amplitude <oo|vv>");
    global_dpd_->buf4_init(&tau1_AO_bb, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[n,n]"), ID("[o,o]"), ID("[n,n]"), 0,
                           "tau1AO <oo|nn>");
    global_dpd_->buf4_scm(&tau1_AO_bb, 0.0);
    half_transform(&tau1_AO_bb, &lambda, bvir_c_, bvir_c_, nbvirpi_, nbvirpi_, pq_row_start, cd_row_start, true, 1.0,
                   0.0);
    global_dpd_->buf4_close(&lambda);
    global_dpd_->buf4_close(&tau1_AO_bb);

    // Now sort for better memory access patterns
    global_dpd_->buf4_init(&tau1_AO_bb, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[n,n]"), ID("[o,o]"), ID("[n,n]"), 1,
                           "tau1AO <oo|nn>");
    global_dpd_->buf4_sort(&tau1_AO_bb, PSIF_DCT_DPD, rspq, ID("[n,n]"), ID("[o,o]"), "tau1AO <nn|oo>");
    global_dpd_->buf4_close(&tau1_AO_bb);

    // Now reopen the two AO dpd_buf4's
    global_dpd_->buf4_init(&tau1_AO_bb, PSIF_DCT_DPD, 0, ID("[n,n]"), ID("[o,o]"), ID("[n,n]"), ID("[o,o]"), 0,
                           "tau1AO <nn|oo>");
    global_dpd_->buf4_init(&tau2_AO_bb, PSIF_DCT_DPD, 0, ID("[n,n]"), ID("[o,o]"), ID("[n,n]"), ID("[o,o]"), 0,
                           "tau2AO <nn|oo>");
    global_dpd_->buf4_scm(&tau2_AO_bb, 0.0);

    /********** AB ***********/
    global_dpd_->buf4_init(&lambda, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0,
                           "Amplitude <Oo|Vv>");
    global_dpd_->buf4_init(&tau1_AO_ab, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[n,n]"), ID("[O,o]"), ID("[n,n]"), 0,
                           "tau1AO <Oo|nn>");
    global_dpd_->buf4_scm(&tau1_AO_ab, 0.0);
    half_transform(&tau1_AO_ab, &lambda, avir_c_, bvir_c_, navirpi_, nbvirpi_, pq_row_start, Cd_row_start, true, 1.0,
                   0.0);
    global_dpd_->buf4_close(&lambda);
    global_dpd_->buf4_close(&tau1_AO_ab);

    // Now sort for better memory access patterns
    global_dpd_->buf4_init(&tau1_AO_ab, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[n,n]"), ID("[O,o]"), ID("[n,n]"), 0,
                           "tau1AO <Oo|nn>");
    global_dpd_->buf4_sort(&tau1_AO_ab, PSIF_DCT_DPD, rspq, ID("[n,n]"), ID("[O,o]"), "tau1AO <nn|Oo>");
    global_dpd_->buf4_close(&tau1_AO_ab);

    // Reopen the two AO dpd_buf4's
    global_dpd_->buf4_init(&tau1_AO_ab, PSIF_DCT_DPD, 0, ID("[n,n]"), ID("[O,o]"), ID("[n,n]"), ID("[O,o]"), 0,
                           "tau1AO <nn|Oo>");
    global_dpd_->buf4_init(&tau2_AO_ab, PSIF_DCT_DPD, 0, ID("[n,n]"), ID("[O,o]"), ID("[n,n]"), ID("[O,o]"), 0,
                           "tau2AO <nn|Oo>");
    global_dpd_->buf4_scm(&tau2_AO_ab, 0.0);

    // Prepare to contract TEI against the 1RDM.

    // Write SO basis gamma to disk as s1(temp)
    global_dpd_->file2_init(&s_aa_1, PSIF_DCT_DPD, 0, ID('n'), ID('n'), "s1(temp)A <n|n>");
    global_dpd_->file2_init(&s_bb_1, PSIF_DCT_DPD, 0, ID('n'), ID('n'), "s1(temp)B <n|n>");
    auto gamma_a = tau_so_a_->clone();
    gamma_a->add(kappa_so_a_);
    auto gamma_b = tau_so_b_->clone();
    gamma_b->add(kappa_so_b_);
    gamma_a->write_to_dpdfile2(&s_aa_1);
    gamma_b->write_to_dpdfile2(&s_bb_1);
    global_dpd_->file2_close(&s_aa_1);
    global_dpd_->file2_close(&s_bb_1);

    // Reopen the arrays
    global_dpd_->file2_init(&s_aa_1, PSIF_DCT_DPD, 0, ID('n'), ID('n'), "s1(temp)A <n|n>");
    global_dpd_->file2_init(&s_bb_1, PSIF_DCT_DPD, 0, ID('n'), ID('n'), "s1(temp)B <n|n>");

    // This is where GGamma contribution will be placed
    global_dpd_->file2_init(&s_aa_2, PSIF_DCT_DPD, 0, ID('n'), ID('n'), "s2(temp)A <n|n>");
    global_dpd_->file2_init(&s_bb_2, PSIF_DCT_DPD, 0, ID('n'), ID('n'), "s2(temp)B <n|n>");
    global_dpd_->file2_scm(&s_aa_2, 0.0);
    global_dpd_->file2_scm(&s_bb_2, 0.0);

    // Now put stuff in memory
    global_dpd_->file2_mat_init(&s_aa_1);
    global_dpd_->file2_mat_init(&s_aa_2);
    global_dpd_->file2_mat_rd(&s_aa_1);
    global_dpd_->file2_mat_rd(&s_aa_2);
    global_dpd_->file2_mat_init(&s_bb_1);
    global_dpd_->file2_mat_init(&s_bb_2);
    global_dpd_->file2_mat_rd(&s_bb_1);
    global_dpd_->file2_mat_rd(&s_bb_2);

    for (int h = 0; h < nirrep_; ++h) {
        global_dpd_->buf4_mat_irrep_init(&tau1_AO_aa, h);
        global_dpd_->buf4_mat_irrep_rd(&tau1_AO_aa, h);
        global_dpd_->buf4_mat_irrep_init(&tau2_AO_aa, h);

        global_dpd_->buf4_mat_irrep_init(&tau1_AO_bb, h);
        global_dpd_->buf4_mat_irrep_rd(&tau1_AO_bb, h);
        global_dpd_->buf4_mat_irrep_init(&tau2_AO_bb, h);

        global_dpd_->buf4_mat_irrep_init(&tau1_AO_ab, h);
        global_dpd_->buf4_mat_irrep_rd(&tau1_AO_ab, h);
        global_dpd_->buf4_mat_irrep_init(&tau2_AO_ab, h);
    }

    bool lastBuffer;
    do {
        lastBuffer = iwl->last_buffer();
        for (int index = 0; index < iwl->buffer_count(); ++index) {
            labelIndex = 4 * index;
            p = std::abs((int)lblptr[labelIndex++]);
            q = (int)lblptr[labelIndex++];
            r = (int)lblptr[labelIndex++];
            s = (int)lblptr[labelIndex++];
            value = (double)valptr[index];
            AO_contribute(&tau1_AO_aa, &tau2_AO_aa, p, q, r, s, value, &s_aa_1, &s_bb_1, &s_aa_2);
            AO_contribute(&tau1_AO_bb, &tau2_AO_bb, p, q, r, s, value, &s_bb_1, &s_aa_1, &s_bb_2);
            AO_contribute(&tau1_AO_ab, &tau2_AO_ab, p, q, r, s, value);
            ++counter;

        } /* end loop through current buffer */
        if (!lastBuffer) iwl->fetch();
    } while (!lastBuffer);
    iwl->set_keep_flag(true);
    delete iwl;
    if (print_ > 1) {
        outfile->Printf("Processed %d SO integrals each for AA, BB, and AB\n", counter);
    }
    for (int h = 0; h < nirrep_; ++h) {
        global_dpd_->buf4_mat_irrep_wrt(&tau2_AO_aa, h);
        global_dpd_->buf4_mat_irrep_close(&tau1_AO_aa, h);
        global_dpd_->buf4_mat_irrep_close(&tau2_AO_aa, h);

        global_dpd_->buf4_mat_irrep_wrt(&tau2_AO_bb, h);
        global_dpd_->buf4_mat_irrep_close(&tau1_AO_bb, h);
        global_dpd_->buf4_mat_irrep_close(&tau2_AO_bb, h);

        global_dpd_->buf4_mat_irrep_wrt(&tau2_AO_ab, h);
        global_dpd_->buf4_mat_irrep_close(&tau1_AO_ab, h);
        global_dpd_->buf4_mat_irrep_close(&tau2_AO_ab, h);
    }

    global_dpd_->buf4_close(&tau1_AO_aa);
    global_dpd_->buf4_close(&tau1_AO_bb);
    global_dpd_->buf4_close(&tau1_AO_ab);
    global_dpd_->buf4_close(&tau2_AO_aa);
    global_dpd_->buf4_close(&tau2_AO_bb);
    global_dpd_->buf4_close(&tau2_AO_ab);
    global_dpd_->file2_mat_wrt(&s_aa_2);
    global_dpd_->file2_mat_wrt(&s_bb_2);
    global_dpd_->file2_mat_close(&s_aa_2);
    global_dpd_->file2_mat_close(&s_bb_2);

    /********** AA ***********/
    global_dpd_->buf4_init(&tau2_AO_aa, PSIF_DCT_DPD, 0, ID("[n,n]"), ID("[O,O]"), ID("[n,n]"), ID("[O,O]"), 0,
                           "tau2AO <nn|OO>");
    global_dpd_->buf4_sort(&tau2_AO_aa, PSIF_DCT_DPD, rspq, ID("[O,O]"), ID("[n,n]"), "tau2AO <OO|nn>");
    global_dpd_->buf4_close(&tau2_AO_aa);
    global_dpd_->buf4_init(&tau2_AO_aa, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[n,n]"), ID("[O,O]"), ID("[n,n]"), 0,
                           "tau2AO <OO|nn>");
    global_dpd_->buf4_init(&tau_temp, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                           "tau(temp) <OO|VV>");
    global_dpd_->buf4_scm(&tau_temp, 0.0);
    half_transform(&tau2_AO_aa, &tau_temp, avir_c_, avir_c_, navirpi_, navirpi_, pq_row_start, CD_row_start, false, 0.5,
                   0.0);
    global_dpd_->buf4_close(&tau2_AO_aa);
    global_dpd_->buf4_close(&tau_temp);

    /********** BB ***********/
    global_dpd_->buf4_init(&tau2_AO_bb, PSIF_DCT_DPD, 0, ID("[n,n]"), ID("[o,o]"), ID("[n,n]"), ID("[o,o]"), 0,
                           "tau2AO <nn|oo>");
    global_dpd_->buf4_sort(&tau2_AO_bb, PSIF_DCT_DPD, rspq, ID("[o,o]"), ID("[n,n]"), "tau2AO <oo|nn>");
    global_dpd_->buf4_close(&tau2_AO_bb);
    global_dpd_->buf4_init(&tau2_AO_bb, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[n,n]"), ID("[o,o]"), ID("[n,n]"), 0,
                           "tau2AO <oo|nn>");
    global_dpd_->buf4_init(&tau_temp, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o,o]"), ID("[v,v]"), 0,
                           "tau(temp) <oo|vv>");
    global_dpd_->buf4_scm(&tau_temp, 0.0);
    half_transform(&tau2_AO_bb, &tau_temp, bvir_c_, bvir_c_, nbvirpi_, nbvirpi_, pq_row_start, cd_row_start, false, 0.5,
                   0.0);
    global_dpd_->buf4_close(&tau2_AO_bb);
    global_dpd_->buf4_close(&tau_temp);

    /********** AB ***********/
    global_dpd_->buf4_init(&tau2_AO_ab, PSIF_DCT_DPD, 0, ID("[n,n]"), ID("[O,o]"), ID("[n,n]"), ID("[O,o]"), 0,
                           "tau2AO <nn|Oo>");
    global_dpd_->buf4_sort(&tau2_AO_ab, PSIF_DCT_DPD, rspq, ID("[O,o]"), ID("[n,n]"), "tau2AO <Oo|nn>");
    global_dpd_->buf4_close(&tau2_AO_ab);
    global_dpd_->buf4_init(&tau2_AO_ab, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[n,n]"), ID("[O,o]"), ID("[n,n]"), 0,
                           "tau2AO <Oo|nn>");
    global_dpd_->buf4_init(&tau_temp, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0,
                           "tau(temp) <Oo|Vv>");
    global_dpd_->buf4_scm(&tau_temp, 0.0);
    half_transform(&tau2_AO_ab, &tau_temp, avir_c_, bvir_c_, navirpi_, nbvirpi_, pq_row_start, Cd_row_start, false, 1.0,
                   0.0);
    global_dpd_->buf4_close(&tau2_AO_ab);
    global_dpd_->buf4_close(&tau_temp);

    // Transform the GGamma contribution from SO to MO basis for the future use in the Fock operator
    // Alpha occupied
    global_dpd_->file2_init(&tau, PSIF_DCT_DPD, 0, ID('O'), ID('O'), "GGamma <O|O>");
    global_dpd_->file2_scm(&tau, 0.0);
    file2_transform(&s_aa_2, &tau, aocc_c_, false);
    global_dpd_->file2_close(&tau);

    // Alpha virtual
    global_dpd_->file2_init(&tau, PSIF_DCT_DPD, 0, ID('V'), ID('V'), "GGamma <V|V>");
    global_dpd_->file2_scm(&tau, 0.0);
    file2_transform(&s_aa_2, &tau, avir_c_, false);
    global_dpd_->file2_close(&tau);

    // Beta occupied
    global_dpd_->file2_init(&tau, PSIF_DCT_DPD, 0, ID('o'), ID('o'), "GGamma <o|o>");
    global_dpd_->file2_scm(&tau, 0.0);
    file2_transform(&s_bb_2, &tau, bocc_c_, false);
    global_dpd_->file2_close(&tau);

    // Beta virtual
    global_dpd_->file2_init(&tau, PSIF_DCT_DPD, 0, ID('v'), ID('v'), "GGamma <v|v>");
    global_dpd_->file2_scm(&tau, 0.0);
    file2_transform(&s_bb_2, &tau, bvir_c_, false);
    global_dpd_->file2_close(&tau);

    global_dpd_->file2_close(&s_aa_1);
    global_dpd_->file2_close(&s_aa_2);
    global_dpd_->file2_close(&s_bb_1);
    global_dpd_->file2_close(&s_bb_2);

    free_int_matrix(pq_row_start);
    free_int_matrix(CD_row_start);
    free_int_matrix(cd_row_start);
    free_int_matrix(Cd_row_start);

    dct_timer_off("DCTSolver::build_tensors");
}

}  // namespace dct
}  // namespace psi
