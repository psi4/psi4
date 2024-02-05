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

#include "dct.h"
#include "psi4/psifiles.h"

#include "psi4/libiwl/iwl.hpp"
#include "psi4/libdpd/dpd.h"
#include "psi4/libqt/qt.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/wavefunction.h"
#include "psi4/libmints/mintshelper.h"
#include "psi4/libtrans/integraltransform.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/liboptions/liboptions.h"

#include <map>
#include <cmath>

namespace psi {
namespace dct {

/**
 * Reads the orbitals and related quantities from the reference wavefunction
 * and reads the one-electron integrals from PSIO.
 * for RHF reference.
 */
void DCTSolver::initialize_orbitals_from_reference_R() {
    dct_timer_on("DCTSolver::rhf_guess");

    epsilon_a_->copy(*reference_wavefunction_->epsilon_a());
    epsilon_b_->copy(*epsilon_a_);
    Ca_->copy(reference_wavefunction_->Ca());
    Cb_->copy(Ca_);
    moFa_->copy(reference_wavefunction_->Fa());
    moFa_->transform(Ca_);
    moFb_->copy(moFa_);
    update_scf_density_RHF();

    dct_timer_off("DCTSolver::rhf_guess");
}

/**
 * Uses the MO coefficients to form the SCF density matrices in the SO basis.
 * @param Whether to damp the update or not
 * @return RMS density change
 */
double DCTSolver::update_scf_density_RHF(bool damp) {
    dct_timer_on("DCTSolver::update_rhf_density");

    double dampingFactor = options_.get_double("DAMPING_PERCENTAGE");  // The default DAMPING_PERCENTAGE is 0.0
    double newFraction = damp ? 1.0 : 1.0 - dampingFactor / 100.0;
    size_t nElements = 0;
    double sumOfSquares = 0.0;
    Matrix old(kappa_so_a_);  // Zero matrix
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

    kappa_so_b_->copy(kappa_so_a_);

    // We're not converged until the RMS error vector *and* the RMS density
    // changes are below the threshold
    dct_timer_off("DCTSolver::update_rhf_density");

    return sqrt(sumOfSquares / nElements);
}

/**
 * Builds the G matrix for Fock matrix, the external potential (tau contracted with the integrals)
 * and other tensors, if requested, out-of-core using the SO integrals.
 * Also builds the AO basis tensors for the <VV||VV>, <vv||vv>, and <Vv|Vv> terms in the G and X intermediates
 * All quantities are built simultaneously to reduce I/O.
 */
void DCTSolver::process_so_ints_RHF() {
    dct_timer_on("DCTSolver::process_so_ints");

    IWL *iwl = new IWL(psio_.get(), PSIF_SO_TEI, int_tolerance_, 1, 1);

    Label *lblptr = iwl->labels();
    Value *valptr = iwl->values();

    double *Da = init_array(ntriso_);
    double *Ga = init_array(ntriso_);

    auto opdm = kappa_so_a_->clone();
    opdm->add(tau_so_a_);
    int soOffset = 0;
    for (int h = 0; h < nirrep_; ++h) {
        for (int mu = 0; mu < nsopi_[h]; ++mu) {
            for (int nu = 0; nu <= mu; ++nu) {
                int muNu = INDEX((nu + soOffset), (mu + soOffset));
                Da[muNu] = opdm->get(h, mu, nu);
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
    dpdbuf4 tau1_AO_ab, tau2_AO_ab;

    bool buildTensors = (options_.get_str("AO_BASIS") == "DISK");

    if (buildTensors) {
        counter = 0;

        // Build the offset arrays needed for the DGEMM in half_transform
        pq_row_start = init_int_matrix(nirrep_, nirrep_);
        Cd_row_start = init_int_matrix(nirrep_, nirrep_);
        for (h = 0; h < nirrep_; ++h) {
            for (Gc = 0, offset = 0; Gc < nirrep_; ++Gc) {
                Gd = Gc ^ h;
                pq_row_start[h][Gc] = offset;
                offset += nsopi_[Gc] * nsopi_[Gd];
            }
            for (Gc = 0, offset = 0; Gc < nirrep_; ++Gc) {
                Gd = Gc ^ h;
                Cd_row_start[h][Gc] = offset;
                offset += navirpi_[Gc] * nbvirpi_[Gd];
            }
        }

        dpd_set_default(_ints->get_dpd_id());

        /********** AB ***********/
        global_dpd_->buf4_init(&lambda, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                               "Amplitude SF <OO|VV>");  // Amplitude <Oo|Vv>
        global_dpd_->buf4_init(&tau1_AO_ab, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[n,n]"), ID("[O,O]"), ID("[n,n]"), 0,
                               "tau1AO <OO|nn>");  // tau1AO <Oo|nn>
        global_dpd_->buf4_scm(&tau1_AO_ab, 0.0);
        half_transform(&tau1_AO_ab, &lambda, avir_c_, avir_c_, navirpi_, navirpi_, pq_row_start, Cd_row_start, true,
                       1.0, 0.0);
        global_dpd_->buf4_close(&lambda);
        global_dpd_->buf4_close(&tau1_AO_ab);

        // Now sort for better memory access patterns
        global_dpd_->buf4_init(&tau1_AO_ab, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[n,n]"), ID("[O,O]"), ID("[n,n]"), 0,
                               "tau1AO <OO|nn>");  // tau1AO <Oo|nn>
        global_dpd_->buf4_sort(&tau1_AO_ab, PSIF_DCT_DPD, rspq, ID("[n,n]"), ID("[O,O]"),
                               "tau1AO <nn|OO>");  // tau1AO <nn|Oo>
        global_dpd_->buf4_close(&tau1_AO_ab);

        // Reopen the two AO dpd_buf4's
        global_dpd_->buf4_init(&tau1_AO_ab, PSIF_DCT_DPD, 0, ID("[n,n]"), ID("[O,O]"), ID("[n,n]"), ID("[O,O]"), 0,
                               "tau1AO <nn|OO>");  // tau1AO <nn|Oo>
        global_dpd_->buf4_init(&tau2_AO_ab, PSIF_DCT_DPD, 0, ID("[n,n]"), ID("[O,O]"), ID("[n,n]"), ID("[O,O]"), 0,
                               "tau2AO <nn|OO>");  // tau2AO <nn|Oo>
        global_dpd_->buf4_scm(&tau2_AO_ab, 0.0);

        // Now put stuff in memory
        for (int h = 0; h < nirrep_; ++h) {
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
            Ga[rsArr] += 2 * Da[pqArr] * value;
            if (q >= r) {
                Ga[qrArr] -= Da[psArr] * value;
            }

            if (p != q && r != s && pqArr != rsArr) {
                /* (pq|sr) */
                if (s >= r) {
                    Ga[srArr] += 2 * Da[pqArr] * value;
                }
                if (q >= s) {
                    Ga[qsArr] -= Da[prArr] * value;
                }

                /* (qp|rs) */
                if (r >= s) {
                    Ga[rsArr] += 2 * Da[qpArr] * value;
                }
                if (p >= r) {
                    Ga[prArr] -= Da[qsArr] * value;
                }

                /* (qp|sr) */
                if (s >= r) {
                    Ga[srArr] += 2 * Da[qpArr] * value;
                }
                if (p >= s) {
                    Ga[psArr] -= Da[qrArr] * value;
                }

                /* (rs|pq) */
                if (p >= q) {
                    Ga[pqArr] += 2 * Da[rsArr] * value;
                }
                if (s >= p) {
                    Ga[spArr] -= Da[rqArr] * value;
                }

                /* (sr|pq) */
                if (p >= q) {
                    Ga[pqArr] += 2 * Da[srArr] * value;
                }
                if (r >= p) {
                    Ga[rpArr] -= Da[sqArr] * value;
                }

                /* (rs|qp) */
                if (q >= p) {
                    Ga[qpArr] += 2 * Da[rsArr] * value;
                }
                if (s >= q) {
                    Ga[sqArr] -= Da[rpArr] * value;
                }

                /* (sr|qp) */
                if (q >= p) {
                    Ga[qpArr] += 2 * Da[srArr] * value;
                }
                if (r >= q) {
                    Ga[rqArr] -= Da[spArr] * value;
                }
            } else if (p != q && r != s && pqArr == rsArr) {
                /* (pq|sr) */
                if (s >= r) {
                    Ga[srArr] += 2 * Da[pqArr] * value;
                }
                if (q >= s) {
                    Ga[qsArr] -= Da[prArr] * value;
                }
                /* (qp|rs) */
                if (r >= s) {
                    Ga[rsArr] += 2 * Da[qpArr] * value;
                }
                if (p >= r) {
                    Ga[prArr] -= Da[qsArr] * value;
                }

                /* (qp|sr) */
                if (s >= r) {
                    Ga[srArr] += 2 * Da[qpArr] * value;
                }
                if (p >= s) {
                    Ga[psArr] -= Da[qrArr] * value;
                }
            } else if (p != q && r == s) {
                /* (qp|rs) */
                if (r >= s) {
                    Ga[rsArr] += 2 * Da[qpArr] * value;
                }
                if (p >= r) {
                    Ga[prArr] -= Da[qsArr] * value;
                }

                /* (rs|pq) */
                if (p >= q) {
                    Ga[pqArr] += 2 * Da[rsArr] * value;
                }
                if (s >= p) {
                    Ga[spArr] -= Da[rqArr] * value;
                }

                /* (rs|qp) */
                if (q >= p) {
                    Ga[qpArr] += 2 * Da[rsArr] * value;
                }
                if (s >= q) {
                    Ga[sqArr] -= Da[rpArr] * value;
                }
            } else if (p == q && r != s) {
                /* (pq|sr) */
                if (s >= r) {
                    Ga[srArr] += 2 * Da[pqArr] * value;
                }
                if (q >= s) {
                    Ga[qsArr] -= Da[prArr] * value;
                }

                /* (rs|pq) */
                if (p >= q) {
                    Ga[pqArr] += 2 * Da[rsArr] * value;
                }
                if (s >= p) {
                    Ga[spArr] -= Da[rqArr] * value;
                }

                /* (sr|pq) */
                if (p >= q) {
                    Ga[pqArr] += 2 * Da[srArr] * value;
                }
                if (r >= p) {
                    Ga[rpArr] -= Da[sqArr] * value;
                }
            } else if (p == q && r == s && pqArr != rsArr) {
                /* (rs|pq) */
                if (p >= q) {
                    Ga[pqArr] += 2 * Da[rsArr] * value;
                }
                if (s >= p) {
                    Ga[spArr] -= Da[rqArr] * value;
                }
            }
        } /* end loop through current buffer */
        if (!lastBuffer) iwl->fetch();
    } while (!lastBuffer);
    iwl->set_keep_flag(true);
    delete iwl;
    if (buildTensors) {
        if (print_ > 1) {
            outfile->Printf("Processed %d SO integrals each for AB\n", counter);
        }
        for (int h = 0; h < nirrep_; ++h) {
            global_dpd_->buf4_mat_irrep_wrt(&tau2_AO_ab, h);
            global_dpd_->buf4_mat_irrep_close(&tau1_AO_ab, h);
            global_dpd_->buf4_mat_irrep_close(&tau2_AO_ab, h);
        }

        global_dpd_->buf4_close(&tau1_AO_ab);
        global_dpd_->buf4_close(&tau2_AO_ab);

        /********** AB ***********/
        global_dpd_->buf4_init(&tau2_AO_ab, PSIF_DCT_DPD, 0, ID("[n,n]"), ID("[O,O]"), ID("[n,n]"), ID("[O,O]"), 0,
                               "tau2AO <nn|OO>");  // tau2AO <nn|Oo>
        global_dpd_->buf4_sort(&tau2_AO_ab, PSIF_DCT_DPD, rspq, ID("[O,O]"), ID("[n,n]"),
                               "tau2AO <OO|nn>");  // tau2AO <Oo|nn>
        global_dpd_->buf4_close(&tau2_AO_ab);
        global_dpd_->buf4_init(&tau2_AO_ab, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[n,n]"), ID("[O,O]"), ID("[n,n]"), 0,
                               "tau2AO <OO|nn>");  // tau2AO <Oo|nn>
        global_dpd_->buf4_init(&tau_temp, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                               "tau(temp) SF <OO|VV>");  // tau(temp) <Oo|Vv>
        global_dpd_->buf4_scm(&tau_temp, 0.0);
        half_transform(&tau2_AO_ab, &tau_temp, avir_c_, avir_c_, navirpi_, navirpi_, pq_row_start, Cd_row_start, false,
                       1.0, 0.0);
        global_dpd_->buf4_close(&tau2_AO_ab);
        global_dpd_->buf4_close(&tau_temp);

        free_int_matrix(pq_row_start);
        free_int_matrix(Cd_row_start);
    }

    // Build the Fock matrices from the H and G matrices
    soOffset = 0;
    for (int h = 0; h < nirrep_; ++h) {
        for (int mu = 0; mu < nsopi_[h]; ++mu) {
            for (int nu = 0; nu <= mu; ++nu) {
                int muNu = INDEX((nu + soOffset), (mu + soOffset));
                double aVal = Ga[muNu];
                Fa_->add(h, mu, nu, aVal);
                if (mu != nu) {
                    Fa_->add(h, nu, mu, aVal);
                }
            }
        }
        soOffset += nsopi_[h];
    }

    Fb_->copy(Fa_);

    free(Da);
    free(Ga);

    dct_timer_off("DCTSolver::process_so_ints");
}

/**
 * Computes the SCF energy from the latest Fock and density matrices.
 * WARNING! This quantity is a misnomer from earlier days of the theory.
 * "SCF" here means "excluding 2RDM cumulant, including 1RDM and 1RDM products."
 */
void DCTSolver::compute_scf_energy_RHF() {
    dct_timer_on("DCTSolver::compute_scf_energy");

    // Escf = eNuc + 0.5 * (H + F) * (kappa + tau)
    scf_energy_ = enuc_;
    scf_energy_ += kappa_so_a_->vector_dot(so_h_);
    scf_energy_ += tau_so_a_->vector_dot(so_h_);

    if (options_.get_str("DCT_TYPE") == "DF" && options_.get_str("AO_BASIS") == "NONE") {
        scf_energy_ += mo_gammaA_.vector_dot(moFa_);
    } else {
        scf_energy_ += kappa_so_a_->vector_dot(Fa_);
        scf_energy_ += tau_so_a_->vector_dot(Fa_);
    }

    dct_timer_off("DCTSolver::compute_scf_energy");
}

double DCTSolver::compute_scf_error_vector_RHF() {
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
    scf_error_b_->copy(scf_error_a_);

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

}  // namespace dct
}  // namespace psi
