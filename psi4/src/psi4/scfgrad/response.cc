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

#include "scf_grad.h"
#include "jk_grad.h"

#include "psi4/libqt/qt.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/libpsio/psio.h"
#include "psi4/psi4-dec.h"
#include "psi4/libmints/dipole.h"
#include "psi4/libmints/cdsalclist.h"
#include "psi4/libfock/v.h"
//#include "psi4/libfock/jk.h"
#include "psi4/libfunctional/superfunctional.h"
#include "psi4/psifiles.h"
#include "psi4/lib3index/dftensor.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/vector.h"
#ifdef USING_ecpint
#include "psi4/libmints/ecpint.h"
#endif
#include "psi4/libmints/integral.h"
#include "psi4/libmints/mintshelper.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/libscf_solver/rhf.h"
#include "psi4/libscf_solver/uhf.h"


#include <algorithm>
#include <mutex>
#include <sstream>

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef USING_BrianQC

#include <use_brian_wrapper.h>
#include <brian_macros.h>
#include <brian_common.h>
#include <brian_cphf.h>

extern void checkBrian();
extern BrianCookie brianCookie;
extern bool brianEnable;
extern bool brianCPHFLeftSideFlag;

#endif

using namespace psi;

namespace psi {
namespace scfgrad {

std::shared_ptr<Matrix> RSCFDeriv::hessian_response() {
    // => Control Parameters <= //

    std::shared_ptr<Vector> eps = epsilon_a_subset("AO", "ALL");
    std::shared_ptr<Vector> eps_occ = epsilon_a_subset("AO", "OCC");
    std::shared_ptr<Vector> eps_vir = epsilon_a_subset("AO", "VIR");
    std::shared_ptr<Matrix> C = Ca_subset("AO", "ALL");
    std::shared_ptr<Matrix> Cocc = Ca_subset("AO", "OCC");
    std::shared_ptr<Matrix> Cvir = Ca_subset("AO", "VIR");
    std::shared_ptr<Matrix> Dt = Da_subset("AO");
    double** Dap = Dt->pointer();

    double** Cp = C->pointer();
    double** Cop = Cocc->pointer();
    double** Cvp = Cvir->pointer();
    double* ep = eps->pointer();
    double* eop = eps_occ->pointer();
    double* evp = eps_vir->pointer();

    // => Sizing <= //

    int natom = molecule_->natom();
    int nso = basisset_->nbf();
    int nocc = eps_occ->dimpi()[0];
    int nvir = eps_vir->dimpi()[0];
    int nmo = C->colspi()[0];

    // => Target <= //

    auto response = std::make_shared<Matrix>("RHF Response", 3 * natom, 3 * natom);

    // => Response Utility File <= //

    psio_->open(PSIF_HESS, PSIO_OPEN_NEW);

    // => Spi <= //
    {
        // Overlap derivatives
        std::shared_ptr<OneBodyAOInt> Sint(integral_->ao_overlap(1));

        auto Smix = std::make_shared<Matrix>("Smix", nso, nocc);
        auto Smiy = std::make_shared<Matrix>("Smiy", nso, nocc);
        auto Smiz = std::make_shared<Matrix>("Smiz", nso, nocc);
        double** Smixp = Smix->pointer();
        double** Smiyp = Smiy->pointer();
        double** Smizp = Smiz->pointer();

        auto Sai = std::make_shared<Matrix>("Sai", nvir, nocc);
        double** Saip = Sai->pointer();
        auto Sij = std::make_shared<Matrix>("Sij", nocc, nocc);
        double** Sijp = Sij->pointer();
        auto Spi = std::make_shared<Matrix>("Spi", nmo, nocc);
        double** Spip = Spi->pointer();

        psio_address next_Sai = PSIO_ZERO;
        psio_address next_Sij = PSIO_ZERO;
        psio_address next_Smi = PSIO_ZERO;
        psio_address next_Spi = PSIO_ZERO;
        for (int A = 0; A < 3 * natom; A++) {
            psio_->write(PSIF_HESS, "Smi^A", (char*)Smixp[0], static_cast<size_t>(nso) * nocc * sizeof(double),
                         next_Smi, &next_Smi);
        }
        for (int A = 0; A < 3 * natom; A++) {
            psio_->write(PSIF_HESS, "Sai^A", (char*)Saip[0], static_cast<size_t>(nvir) * nocc * sizeof(double),
                         next_Sai, &next_Sai);
        }
        for (int A = 0; A < 3 * natom; A++) {
            psio_->write(PSIF_HESS, "Sij^A", (char*)Sijp[0], static_cast<size_t>(nocc) * nocc * sizeof(double),
                         next_Sij, &next_Sij);
        }
        for (int A = 0; A < 3 * natom; A++) {
            psio_->write(PSIF_HESS, "Spi^A", (char*)Spip[0], static_cast<size_t>(nmo) * nocc * sizeof(double), next_Spi,
                         &next_Spi);
        }
        next_Smi = PSIO_ZERO;
        next_Sai = PSIO_ZERO;
        next_Sij = PSIO_ZERO;
        next_Spi = PSIO_ZERO;

#ifdef USING_BrianQC
        brianInt maxSegmentSize;
        brianInt maxSegmentAtomCount;
        brianInt segmentAtomCount;
        brianInt segmentAtomIndexStart;
        
        std::vector<std::shared_ptr<Matrix>> Smnx;
        std::vector<std::shared_ptr<Matrix>> Smny;
        std::vector<std::shared_ptr<Matrix>> Smnz;
        std::vector<double*> Smn;
        
        if (brianEnable) {
            brianCPHFMaxSegmentSize(&brianCookie, &maxSegmentSize);
            maxSegmentAtomCount = maxSegmentSize / 3;
            segmentAtomCount = -1;
            segmentAtomIndexStart = -1;
            
            Smnx.resize(maxSegmentAtomCount);
            Smny.resize(maxSegmentAtomCount);
            Smnz.resize(maxSegmentAtomCount);
            Smn.resize(maxSegmentAtomCount * 3);
            for (int i = 0; i < maxSegmentAtomCount; i++) {
                Smnx[i] = std::make_shared<Matrix>("Smnx", nso, nso);
                Smny[i] = std::make_shared<Matrix>("Smny", nso, nso);
                Smnz[i] = std::make_shared<Matrix>("Smnz", nso, nso);
                Smn[i * 3 + 0] = Smnx[i]->get_pointer();
                Smn[i * 3 + 1] = Smny[i]->get_pointer();
                Smn[i * 3 + 2] = Smnz[i]->get_pointer();
            }
        }
#endif

        for (int A = 0; A < natom; A++) {
#ifdef USING_BrianQC
            if (brianEnable) {
                if (segmentAtomCount < 0 || A < segmentAtomIndexStart || A >= (segmentAtomIndexStart + segmentAtomCount)) {
                    segmentAtomIndexStart = A;
                    segmentAtomCount = (segmentAtomIndexStart + maxSegmentAtomCount > natom) ? (natom - segmentAtomIndexStart) : maxSegmentAtomCount;
                    
                    brianInt integralType = BRIAN_INTEGRAL_TYPE_OVERLAP;
                    brianCPHFBuild1eDeriv(&brianCookie, &integralType, &segmentAtomCount, &segmentAtomIndexStart, Smn.data());
                }
                
                C_DGEMM('N', 'N', nso, nocc, nso, 2.0, Smnx[A - segmentAtomIndexStart]->get_pointer(), nso, Cocc->get_pointer(), nocc, 0.0, Smix->get_pointer(), nocc);
                C_DGEMM('N', 'N', nso, nocc, nso, 2.0, Smny[A - segmentAtomIndexStart]->get_pointer(), nso, Cocc->get_pointer(), nocc, 0.0, Smiy->get_pointer(), nocc);
                C_DGEMM('N', 'N', nso, nocc, nso, 2.0, Smnz[A - segmentAtomIndexStart]->get_pointer(), nso, Cocc->get_pointer(), nocc, 0.0, Smiz->get_pointer(), nocc);
            } else {
#endif
            Smix->zero();
            Smiy->zero();
            Smiz->zero();
            const auto& shell_pairs = Sint->shellpairs();
            size_t n_pairs = shell_pairs.size();

            for (size_t p = 0; p < n_pairs; ++p) {
                auto P = shell_pairs[p].first;
                auto Q = shell_pairs[p].second;
                const auto &shellP = basisset_->shell(P);
                const auto &shellQ = basisset_->shell(Q);
                int aP = shellP.ncenter();
                int aQ = shellQ.ncenter();
                if ((aP != A && aQ != A) || aP == aQ) continue;
                Sint->compute_shell_deriv1(P, Q);
                const auto &buffers = Sint->buffers();
                int nP = shellP.nfunction();
                int nQ = shellQ.nfunction();
                int oP = shellP.function_index();
                int oQ = shellQ.function_index();
                const double* buffer2;
                double scale = P == Q ? 1.0 : 2.0;

                if (aP == A) {
                    // Px
                    buffer2 = buffers[0];
                    for (int p = 0; p < nP; p++) {
                        for (int q = 0; q < nQ; q++) {
                            C_DAXPY(nocc, scale * (*buffer2), Cop[q + oQ], 1, Smixp[p + oP], 1);
                            C_DAXPY(nocc, scale * (*buffer2++), Cop[p + oP], 1, Smixp[q + oQ], 1);
                        }
                    }
                    // Py
                    buffer2 = buffers[1];
                    for (int p = 0; p < nP; p++) {
                        for (int q = 0; q < nQ; q++) {
                            C_DAXPY(nocc, scale * (*buffer2), Cop[q + oQ], 1, Smiyp[p + oP], 1);
                            C_DAXPY(nocc, scale * (*buffer2++), Cop[p + oP], 1, Smiyp[q + oQ], 1);
                        }
                    }
                    // Pz
                    buffer2 = buffers[2];
                    for (int p = 0; p < nP; p++) {
                        for (int q = 0; q < nQ; q++) {
                            C_DAXPY(nocc, scale * (*buffer2), Cop[q + oQ], 1, Smizp[p + oP], 1);
                            C_DAXPY(nocc, scale * (*buffer2++), Cop[p + oP], 1, Smizp[q + oQ], 1);
                        }
                    }
                }
                if (aQ == A) {
                    // Qx
                    buffer2 = buffers[3];
                    for (int p = 0; p < nP; p++) {
                        for (int q = 0; q < nQ; q++) {
                            C_DAXPY(nocc, scale * (*buffer2), Cop[q + oQ], 1, Smixp[p + oP], 1);
                            C_DAXPY(nocc, scale * (*buffer2++), Cop[p + oP], 1, Smixp[q + oQ], 1);
                        }
                    }
                    // Qy
                    buffer2 = buffers[4];
                    for (int p = 0; p < nP; p++) {
                        for (int q = 0; q < nQ; q++) {
                            C_DAXPY(nocc, scale * (*buffer2), Cop[q + oQ], 1, Smiyp[p + oP], 1);
                            C_DAXPY(nocc, scale * (*buffer2++), Cop[p + oP], 1, Smiyp[q + oQ], 1);
                        }
                    }
                    // Qz
                    buffer2 = buffers[5];
                    for (int p = 0; p < nP; p++) {
                        for (int q = 0; q < nQ; q++) {
                            C_DAXPY(nocc, scale * (*buffer2), Cop[q + oQ], 1, Smizp[p + oP], 1);
                            C_DAXPY(nocc, scale * (*buffer2++), Cop[p + oP], 1, Smizp[q + oQ], 1);
                        }
                    }
                }
            }
#ifdef USING_BrianQC
            }
#endif
            // Smi_x
            psio_->write(PSIF_HESS, "Smi^A", (char*)Smixp[0], static_cast<size_t>(nso) * nocc * sizeof(double),
                         next_Smi, &next_Smi);
            // Smi_y
            psio_->write(PSIF_HESS, "Smi^A", (char*)Smiyp[0], static_cast<size_t>(nso) * nocc * sizeof(double),
                         next_Smi, &next_Smi);
            // Smi_z
            psio_->write(PSIF_HESS, "Smi^A", (char*)Smizp[0], static_cast<size_t>(nso) * nocc * sizeof(double),
                         next_Smi, &next_Smi);

            // Sai_x
            C_DGEMM('T', 'N', nvir, nocc, nso, 0.5, Cvp[0], nvir, Smixp[0], nocc, 0.0, Saip[0], nocc);
            psio_->write(PSIF_HESS, "Sai^A", (char*)Saip[0], static_cast<size_t>(nvir) * nocc * sizeof(double),
                         next_Sai, &next_Sai);
            // Sai_y
            C_DGEMM('T', 'N', nvir, nocc, nso, 0.5, Cvp[0], nvir, Smiyp[0], nocc, 0.0, Saip[0], nocc);
            psio_->write(PSIF_HESS, "Sai^A", (char*)Saip[0], static_cast<size_t>(nvir) * nocc * sizeof(double),
                         next_Sai, &next_Sai);
            // Sai_z
            C_DGEMM('T', 'N', nvir, nocc, nso, 0.5, Cvp[0], nvir, Smizp[0], nocc, 0.0, Saip[0], nocc);
            psio_->write(PSIF_HESS, "Sai^A", (char*)Saip[0], static_cast<size_t>(nvir) * nocc * sizeof(double),
                         next_Sai, &next_Sai);

            // Sij_x
            C_DGEMM('T', 'N', nocc, nocc, nso, 0.5, Cop[0], nocc, Smixp[0], nocc, 0.0, Sijp[0], nocc);
            psio_->write(PSIF_HESS, "Sij^A", (char*)Sijp[0], static_cast<size_t>(nocc) * nocc * sizeof(double),
                         next_Sij, &next_Sij);
            // Sij_y
            C_DGEMM('T', 'N', nocc, nocc, nso, 0.5, Cop[0], nocc, Smiyp[0], nocc, 0.0, Sijp[0], nocc);
            psio_->write(PSIF_HESS, "Sij^A", (char*)Sijp[0], static_cast<size_t>(nocc) * nocc * sizeof(double),
                         next_Sij, &next_Sij);
            // Sij_z
            C_DGEMM('T', 'N', nocc, nocc, nso, 0.5, Cop[0], nocc, Smizp[0], nocc, 0.0, Sijp[0], nocc);
            psio_->write(PSIF_HESS, "Sij^A", (char*)Sijp[0], static_cast<size_t>(nocc) * nocc * sizeof(double),
                         next_Sij, &next_Sij);

            // Spi_x
            C_DGEMM('T', 'N', nmo, nocc, nso, 0.5, Cp[0], nmo, Smixp[0], nocc, 0.0, Spip[0], nocc);
            psio_->write(PSIF_HESS, "Spi^A", (char*)Spip[0], static_cast<size_t>(nmo) * nocc * sizeof(double), next_Spi,
                         &next_Spi);
            // Spi_y
            C_DGEMM('T', 'N', nmo, nocc, nso, 0.5, Cp[0], nmo, Smiyp[0], nocc, 0.0, Spip[0], nocc);
            psio_->write(PSIF_HESS, "Spi^A", (char*)Spip[0], static_cast<size_t>(nmo) * nocc * sizeof(double), next_Spi,
                         &next_Spi);
            // Spi_z
            C_DGEMM('T', 'N', nmo, nocc, nso, 0.5, Cp[0], nmo, Smizp[0], nocc, 0.0, Spip[0], nocc);
            psio_->write(PSIF_HESS, "Spi^A", (char*)Spip[0], static_cast<size_t>(nmo) * nocc * sizeof(double), next_Spi,
                         &next_Spi);
        }
    }

    // => Tpi <= //
    {
        // Kinetic derivatives
        std::shared_ptr<OneBodyAOInt> Tint(integral_->ao_kinetic(1));

        auto Tmix = std::make_shared<Matrix>("Tmix", nso, nocc);
        auto Tmiy = std::make_shared<Matrix>("Tmiy", nso, nocc);
        auto Tmiz = std::make_shared<Matrix>("Tmiz", nso, nocc);
        double** Tmixp = Tmix->pointer();
        double** Tmiyp = Tmiy->pointer();
        double** Tmizp = Tmiz->pointer();

        auto Tpi = std::make_shared<Matrix>("Tpi", nmo, nocc);
        double** Tpip = Tpi->pointer();
        psio_address next_Tpi = PSIO_ZERO;

#ifdef USING_BrianQC
        brianInt maxSegmentSize;
        brianInt maxSegmentAtomCount;
        brianInt segmentAtomCount;
        brianInt segmentAtomIndexStart;
        
        std::vector<std::shared_ptr<Matrix>> Tmnx;
        std::vector<std::shared_ptr<Matrix>> Tmny;
        std::vector<std::shared_ptr<Matrix>> Tmnz;
        std::vector<double*> Tmn;
        
        if (brianEnable) {
            brianCPHFMaxSegmentSize(&brianCookie, &maxSegmentSize);
            maxSegmentAtomCount = maxSegmentSize / 3;
            segmentAtomCount = -1;
            segmentAtomIndexStart = -1;
            
            Tmnx.resize(maxSegmentAtomCount);
            Tmny.resize(maxSegmentAtomCount);
            Tmnz.resize(maxSegmentAtomCount);
            Tmn.resize(maxSegmentAtomCount * 3);
            for (int i = 0; i < maxSegmentAtomCount; i++) {
                Tmnx[i] = std::make_shared<Matrix>("Tmnx", nso, nso);
                Tmny[i] = std::make_shared<Matrix>("Tmny", nso, nso);
                Tmnz[i] = std::make_shared<Matrix>("Tmnz", nso, nso);
                Tmn[i * 3 + 0] = Tmnx[i]->get_pointer();
                Tmn[i * 3 + 1] = Tmny[i]->get_pointer();
                Tmn[i * 3 + 2] = Tmnz[i]->get_pointer();
            }
        }
#endif

        for (int A = 0; A < natom; A++) {
#ifdef USING_BrianQC
            if (brianEnable) {
                if (segmentAtomCount < 0 || A < segmentAtomIndexStart || A >= (segmentAtomIndexStart + segmentAtomCount)) {
                    segmentAtomIndexStart = A;
                    segmentAtomCount = (segmentAtomIndexStart + maxSegmentAtomCount > natom) ? (natom - segmentAtomIndexStart) : maxSegmentAtomCount;
                    
                    brianInt integralType = BRIAN_INTEGRAL_TYPE_KINETIC;
                    brianCPHFBuild1eDeriv(&brianCookie, &integralType, &segmentAtomCount, &segmentAtomIndexStart, Tmn.data());
                }
                
                C_DGEMM('N', 'N', nso, nocc, nso, 2.0, Tmnx[A - segmentAtomIndexStart]->get_pointer(), nso, Cocc->get_pointer(), nocc, 0.0, Tmix->get_pointer(), nocc);
                C_DGEMM('N', 'N', nso, nocc, nso, 2.0, Tmny[A - segmentAtomIndexStart]->get_pointer(), nso, Cocc->get_pointer(), nocc, 0.0, Tmiy->get_pointer(), nocc);
                C_DGEMM('N', 'N', nso, nocc, nso, 2.0, Tmnz[A - segmentAtomIndexStart]->get_pointer(), nso, Cocc->get_pointer(), nocc, 0.0, Tmiz->get_pointer(), nocc);
            } else {
#endif
            Tmix->zero();
            Tmiy->zero();
            Tmiz->zero();
            const auto& shell_pairs = Tint->shellpairs();
            size_t n_pairs = shell_pairs.size();

            for (size_t p = 0; p < n_pairs; ++p) {
                auto P = shell_pairs[p].first;
                auto Q = shell_pairs[p].second;
                const auto & shellP = basisset_->shell(P);
                const auto & shellQ = basisset_->shell(Q);
                int aP = shellP.ncenter();
                int aQ = shellQ.ncenter();
                if ((aP != A && aQ != A) || aP == aQ) continue;
                Tint->compute_shell_deriv1(P, Q);
                const auto &buffers = Tint->buffers();
                int nP = shellP.nfunction();
                int nQ = shellQ.nfunction();
                int oP = shellP.function_index();
                int oQ = shellQ.function_index();
                const double* buffer2;

                double scale = P == Q ? 1.0 : 2.0;
                if (aP == A) {
                    // Px
                    buffer2 = buffers[0];
                    for (int p = 0; p < nP; p++) {
                        for (int q = 0; q < nQ; q++) {
                            C_DAXPY(nocc, scale * (*buffer2), Cop[q + oQ], 1, Tmixp[p + oP], 1);
                            C_DAXPY(nocc, scale * (*buffer2++), Cop[p + oP], 1, Tmixp[q + oQ], 1);
                        }
                    }
                    // Py
                    buffer2 = buffers[1];
                    for (int p = 0; p < nP; p++) {
                        for (int q = 0; q < nQ; q++) {
                            C_DAXPY(nocc, scale * (*buffer2), Cop[q + oQ], 1, Tmiyp[p + oP], 1);
                            C_DAXPY(nocc, scale * (*buffer2++), Cop[p + oP], 1, Tmiyp[q + oQ], 1);
                        }
                    }
                    // Pz
                    buffer2 = buffers[2];
                    for (int p = 0; p < nP; p++) {
                        for (int q = 0; q < nQ; q++) {
                            C_DAXPY(nocc, scale * (*buffer2), Cop[q + oQ], 1, Tmizp[p + oP], 1);
                            C_DAXPY(nocc, scale * (*buffer2++), Cop[p + oP], 1, Tmizp[q + oQ], 1);
                        }
                    }
                }
                if (aQ == A) {
                    // Qx
                    buffer2 = buffers[3];
                    for (int p = 0; p < nP; p++) {
                        for (int q = 0; q < nQ; q++) {
                            C_DAXPY(nocc, scale * (*buffer2), Cop[q + oQ], 1, Tmixp[p + oP], 1);
                            C_DAXPY(nocc, scale * (*buffer2++), Cop[p + oP], 1, Tmixp[q + oQ], 1);
                        }
                    }
                    // Qy
                    buffer2 = buffers[4];
                    for (int p = 0; p < nP; p++) {
                        for (int q = 0; q < nQ; q++) {
                            C_DAXPY(nocc, scale * (*buffer2), Cop[q + oQ], 1, Tmiyp[p + oP], 1);
                            C_DAXPY(nocc, scale * (*buffer2++), Cop[p + oP], 1, Tmiyp[q + oQ], 1);
                        }
                    }
                    // Qz
                    buffer2 = buffers[5];
                    for (int p = 0; p < nP; p++) {
                        for (int q = 0; q < nQ; q++) {
                            C_DAXPY(nocc, scale * (*buffer2), Cop[q + oQ], 1, Tmizp[p + oP], 1);
                            C_DAXPY(nocc, scale * (*buffer2++), Cop[p + oP], 1, Tmizp[q + oQ], 1);
                        }
                    }
                }
            }
#ifdef USING_BrianQC
            }
#endif
            // Tpi_x
            C_DGEMM('T', 'N', nmo, nocc, nso, 0.5, Cp[0], nmo, Tmixp[0], nocc, 0.0, Tpip[0], nocc);
            psio_->write(PSIF_HESS, "Tpi^A", (char*)Tpip[0], static_cast<size_t>(nmo) * nocc * sizeof(double), next_Tpi,
                         &next_Tpi);
            // Tpi_y
            C_DGEMM('T', 'N', nmo, nocc, nso, 0.5, Cp[0], nmo, Tmiyp[0], nocc, 0.0, Tpip[0], nocc);
            psio_->write(PSIF_HESS, "Tpi^A", (char*)Tpip[0], static_cast<size_t>(nmo) * nocc * sizeof(double), next_Tpi,
                         &next_Tpi);
            // Tpi_z
            C_DGEMM('T', 'N', nmo, nocc, nso, 0.5, Cp[0], nmo, Tmizp[0], nocc, 0.0, Tpip[0], nocc);
            psio_->write(PSIF_HESS, "Tpi^A", (char*)Tpip[0], static_cast<size_t>(nmo) * nocc * sizeof(double), next_Tpi,
                         &next_Tpi);
        }
    }

    if (basisset_->has_ECP()) {
    // => ECP <= //
#ifdef USING_ecpint
    {
        // Effective core potential derivatives
        std::shared_ptr<ECPInt> ecpint(dynamic_cast<ECPInt*>(integral_->ao_ecp(1).release()));

        auto Emix = std::make_shared<Matrix>("Emix", nso, nocc);
        auto Emiy = std::make_shared<Matrix>("Emiy", nso, nocc);
        auto Emiz = std::make_shared<Matrix>("Emiz", nso, nocc);
        double** Emixp = Emix->pointer();
        double** Emiyp = Emiy->pointer();
        double** Emizp = Emiz->pointer();

        // Make a list of all ECP centers
        std::set<int> ecp_centers;
        for (int ecp_shell = 0; ecp_shell < basisset_->n_ecp_shell(); ++ecp_shell){
            const GaussianShell &ecp = basisset_->ecp_shell(ecp_shell);
            ecp_centers.insert(ecp.ncenter());
        }

        auto Epi = std::make_shared<Matrix>("Epi", nmo, nocc);
        double** Epip = Epi->pointer();
        psio_address next_Epi = PSIO_ZERO;


        for (int A = 0; A < natom; A++) {
            Emix->zero();
            Emiy->zero();
            Emiz->zero();
            const auto& shell_pairs = ecpint->shellpairs();
            size_t n_pairs = shell_pairs.size();

            for (size_t p = 0; p < n_pairs; ++p) {
                auto P = shell_pairs[p].first;
                auto Q = shell_pairs[p].second;
                const auto & shellP = basisset_->shell(P);
                const auto & shellQ = basisset_->shell(Q);
                int aP = shellP.ncenter();
                int aQ = shellQ.ncenter();

                // Make a list of all ECP centers and the current basis function pair
                std::set<int> all_centers(ecp_centers.begin(), ecp_centers.end());
                all_centers.insert(aP);
                all_centers.insert(aQ);
                if (all_centers.find(A) == all_centers.end()) continue;

                ecpint->compute_shell_deriv1(P, Q);
                const auto &buffers = ecpint->buffers();
                int nP = shellP.nfunction();
                int nQ = shellQ.nfunction();
                int oP = shellP.function_index();
                int oQ = shellQ.function_index();
                const double* buffer2;

                for (const int center : all_centers) {
                    double scale = P == Q ? 1.0 : 2.0;
                    if (center == A) {
                        // x
                        buffer2 = buffers[3*center+0];
                        for (int p = 0; p < nP; p++) {
                            for (int q = 0; q < nQ; q++) {
                                C_DAXPY(nocc, scale * (*buffer2), Cop[q + oQ], 1, Emixp[p + oP], 1);
                                C_DAXPY(nocc, scale * (*buffer2++), Cop[p + oP], 1, Emixp[q + oQ], 1);
                            }
                        }
                        // y
                        buffer2 = buffers[3*center+1];
                        for (int p = 0; p < nP; p++) {
                            for (int q = 0; q < nQ; q++) {
                                C_DAXPY(nocc, scale * (*buffer2), Cop[q + oQ], 1, Emiyp[p + oP], 1);
                                C_DAXPY(nocc, scale * (*buffer2++), Cop[p + oP], 1, Emiyp[q + oQ], 1);
                            }
                        }
                        // z
                        buffer2 = buffers[3*center+2];
                        for (int p = 0; p < nP; p++) {
                            for (int q = 0; q < nQ; q++) {
                                C_DAXPY(nocc, scale * (*buffer2), Cop[q + oQ], 1, Emizp[p + oP], 1);
                                C_DAXPY(nocc, scale * (*buffer2++), Cop[p + oP], 1, Emizp[q + oQ], 1);
                            }
                        }
                    }
                }
            }
            // Epi_x
            C_DGEMM('T', 'N', nmo, nocc, nso, 0.5, Cp[0], nmo, Emixp[0], nocc, 0.0, Epip[0], nocc);
            psio_->write(PSIF_HESS, "Epi^A", (char*)Epip[0], static_cast<size_t>(nmo) * nocc * sizeof(double), next_Epi,
                         &next_Epi);
            // Epi_y
            C_DGEMM('T', 'N', nmo, nocc, nso, 0.5, Cp[0], nmo, Emiyp[0], nocc, 0.0, Epip[0], nocc);
            psio_->write(PSIF_HESS, "Epi^A", (char*)Epip[0], static_cast<size_t>(nmo) * nocc * sizeof(double), next_Epi,
                         &next_Epi);
            // Epi_z
            C_DGEMM('T', 'N', nmo, nocc, nso, 0.5, Cp[0], nmo, Emizp[0], nocc, 0.0, Epip[0], nocc);
            psio_->write(PSIF_HESS, "Epi^A", (char*)Epip[0], static_cast<size_t>(nmo) * nocc * sizeof(double), next_Epi,
                         &next_Epi);
        }
    }
#endif
    }

    // => Vpi <= //
    {
        // Potential derivatives
        std::shared_ptr<OneBodyAOInt> Vint(integral_->ao_potential(1));

        auto Vmix = std::make_shared<Matrix>("Vmix", nso, nocc);
        auto Vmiy = std::make_shared<Matrix>("Vmiy", nso, nocc);
        auto Vmiz = std::make_shared<Matrix>("Vmiz", nso, nocc);
        double** Vmixp = Vmix->pointer();
        double** Vmiyp = Vmiy->pointer();
        double** Vmizp = Vmiz->pointer();

        auto Vpi = std::make_shared<Matrix>("Vpi", nmo, nocc);
        double** Vpip = Vpi->pointer();
        psio_address next_Vpi = PSIO_ZERO;

#ifdef USING_BrianQC
        brianInt maxSegmentSize;
        brianInt maxSegmentAtomCount;
        brianInt segmentAtomCount;
        brianInt segmentAtomIndexStart;
        
        std::vector<std::shared_ptr<Matrix>> Vmnx;
        std::vector<std::shared_ptr<Matrix>> Vmny;
        std::vector<std::shared_ptr<Matrix>> Vmnz;
        std::vector<double*> Vmn;
        
        if (brianEnable) {
            brianCPHFMaxSegmentSize(&brianCookie, &maxSegmentSize);
            maxSegmentAtomCount = maxSegmentSize / 3;
            segmentAtomCount = -1;
            segmentAtomIndexStart = -1;
            
            Vmnx.resize(maxSegmentAtomCount);
            Vmny.resize(maxSegmentAtomCount);
            Vmnz.resize(maxSegmentAtomCount);
            Vmn.resize(maxSegmentAtomCount * 3);
            for (int i = 0; i < maxSegmentAtomCount; i++) {
                Vmnx[i] = std::make_shared<Matrix>("Vmnx", nso, nso);
                Vmny[i] = std::make_shared<Matrix>("Vmny", nso, nso);
                Vmnz[i] = std::make_shared<Matrix>("Vmnz", nso, nso);
                Vmn[i * 3 + 0] = Vmnx[i]->get_pointer();
                Vmn[i * 3 + 1] = Vmny[i]->get_pointer();
                Vmn[i * 3 + 2] = Vmnz[i]->get_pointer();
            }
        }
#endif

        const auto& shell_pairs = Vint->shellpairs();
        size_t n_pairs = shell_pairs.size();
        for (int A = 0; A < natom; A++) {
#ifdef USING_BrianQC
            if (brianEnable) {
                if (segmentAtomCount < 0 || A < segmentAtomIndexStart || A >= (segmentAtomIndexStart + segmentAtomCount)) {
                    segmentAtomIndexStart = A;
                    segmentAtomCount = (segmentAtomIndexStart + maxSegmentAtomCount > natom) ? (natom - segmentAtomIndexStart) : maxSegmentAtomCount;
                    
                    brianInt integralType = BRIAN_INTEGRAL_TYPE_NUCLEAR;
                    brianCPHFBuild1eDeriv(&brianCookie, &integralType, &segmentAtomCount, &segmentAtomIndexStart, Vmn.data());
                }
                
                C_DGEMM('N', 'N', nso, nocc, nso, 1.0, Vmnx[A - segmentAtomIndexStart]->get_pointer(), nso, Cocc->get_pointer(), nocc, 0.0, Vmix->get_pointer(), nocc);
                C_DGEMM('N', 'N', nso, nocc, nso, 1.0, Vmny[A - segmentAtomIndexStart]->get_pointer(), nso, Cocc->get_pointer(), nocc, 0.0, Vmiy->get_pointer(), nocc);
                C_DGEMM('N', 'N', nso, nocc, nso, 1.0, Vmnz[A - segmentAtomIndexStart]->get_pointer(), nso, Cocc->get_pointer(), nocc, 0.0, Vmiz->get_pointer(), nocc);
            } else {
#endif
            Vmix->zero();
            Vmiy->zero();
            Vmiz->zero();

            for (size_t p = 0; p < n_pairs; ++p) {
                auto P = shell_pairs[p].first;
                auto Q = shell_pairs[p].second;
                const auto & shellP = basisset_->shell(P);
                const auto & shellQ = basisset_->shell(Q);
                int aP = shellP.ncenter();
                int aQ = shellQ.ncenter();
                Vint->compute_shell_deriv1(P, Q);
                const auto &buffers = Vint->buffers();
                int nP = shellP.nfunction();
                int nQ = shellQ.nfunction();
                int oP = shellP.function_index();
                int oQ = shellQ.function_index();

                // buffer ordering is [Px, Py, Pz, Qx, Qy, Qz, A1x, A1y, A1z, A2x, ... ANy, ANz] where
                // the Px is the x derivative of shell P and A1x is the derivative w.r.t. the position
                // of the nuclear charge located on atom1.  There are therefore 6 + 3*natoms buffers.
                const double* buf_x = buffers[3 * A + 6];
                const double* buf_y = buffers[3 * A + 7];
                const double* buf_z = buffers[3 * A + 8];
                double scale = P == Q ? 0.5 : 1.0;

                // Ax
                for (int p = 0; p < nP; p++) {
                    for (int q = 0; q < nQ; q++) {
                        C_DAXPY(nocc, scale * (*buf_x), Cop[q + oQ], 1, Vmixp[p + oP], 1);
                        C_DAXPY(nocc, scale * (*buf_x++), Cop[p + oP], 1, Vmixp[q + oQ], 1);
                    }
                }

                // Ay
                for (int p = 0; p < nP; p++) {
                    for (int q = 0; q < nQ; q++) {
                        C_DAXPY(nocc, scale * (*buf_y), Cop[q + oQ], 1, Vmiyp[p + oP], 1);
                        C_DAXPY(nocc, scale * (*buf_y++), Cop[p + oP], 1, Vmiyp[q + oQ], 1);
                    }
                }

                // Az
                for (int p = 0; p < nP; p++) {
                    for (int q = 0; q < nQ; q++) {
                        C_DAXPY(nocc, scale * (*buf_z), Cop[q + oQ], 1, Vmizp[p + oP], 1);
                        C_DAXPY(nocc, scale * (*buf_z++), Cop[p + oP], 1, Vmizp[q + oQ], 1);
                    }
                }

                // (P| derivatives
                if (aP == A) {
                    const double* buf_x = buffers[0];
                    const double* buf_y = buffers[1];
                    const double* buf_z = buffers[2];
                    // Px
                    for (int p = 0; p < nP; p++) {
                        for (int q = 0; q < nQ; q++) {
                            C_DAXPY(nocc, scale * (*buf_x), Cop[q + oQ], 1, Vmixp[p + oP], 1);
                            C_DAXPY(nocc, scale * (*buf_x++), Cop[p + oP], 1, Vmixp[q + oQ], 1);
                        }
                    }

                    // Py
                    for (int p = 0; p < nP; p++) {
                        for (int q = 0; q < nQ; q++) {
                            C_DAXPY(nocc, scale * (*buf_y), Cop[q + oQ], 1, Vmiyp[p + oP], 1);
                            C_DAXPY(nocc, scale * (*buf_y++), Cop[p + oP], 1, Vmiyp[q + oQ], 1);
                        }
                    }

                    // Pz
                    for (int p = 0; p < nP; p++) {
                        for (int q = 0; q < nQ; q++) {
                            C_DAXPY(nocc, scale * (*buf_z), Cop[q + oQ], 1, Vmizp[p + oP], 1);
                            C_DAXPY(nocc, scale * (*buf_z++), Cop[p + oP], 1, Vmizp[q + oQ], 1);
                        }
                    }
                }

                // |Q) derivatives
                if (aQ == A) {
                    const double* buf_x = buffers[3];
                    const double* buf_y = buffers[4];
                    const double* buf_z = buffers[5];
                    // Qx
                    for (int p = 0; p < nP; p++) {
                        for (int q = 0; q < nQ; q++) {
                            C_DAXPY(nocc, scale * (*buf_x), Cop[q + oQ], 1, Vmixp[p + oP], 1);
                            C_DAXPY(nocc, scale * (*buf_x++), Cop[p + oP], 1, Vmixp[q + oQ], 1);
                        }
                    }

                    // Qy
                    for (int p = 0; p < nP; p++) {
                        for (int q = 0; q < nQ; q++) {
                            C_DAXPY(nocc, scale * (*buf_y), Cop[q + oQ], 1, Vmiyp[p + oP], 1);
                            C_DAXPY(nocc, scale * (*buf_y++), Cop[p + oP], 1, Vmiyp[q + oQ], 1);
                        }
                    }

                    // Qz
                    for (int p = 0; p < nP; p++) {
                        for (int q = 0; q < nQ; q++) {
                            C_DAXPY(nocc, scale * (*buf_z), Cop[q + oQ], 1, Vmizp[p + oP], 1);
                            C_DAXPY(nocc, scale * (*buf_z++), Cop[p + oP], 1, Vmizp[q + oQ], 1);
                        }
                    }
                }
            }
#ifdef USING_BrianQC
            }
#endif
            // Vpi_x
            C_DGEMM('T', 'N', nmo, nocc, nso, 1.0, Cp[0], nmo, Vmixp[0], nocc, 0.0, Vpip[0], nocc);
            psio_->write(PSIF_HESS, "Vpi^A", (char*)Vpip[0], static_cast<size_t>(nmo) * nocc * sizeof(double), next_Vpi,
                         &next_Vpi);
            // Vpi_y
            C_DGEMM('T', 'N', nmo, nocc, nso, 1.0, Cp[0], nmo, Vmiyp[0], nocc, 0.0, Vpip[0], nocc);
            psio_->write(PSIF_HESS, "Vpi^A", (char*)Vpip[0], static_cast<size_t>(nmo) * nocc * sizeof(double), next_Vpi,
                         &next_Vpi);
            // Vpi_z
            C_DGEMM('T', 'N', nmo, nocc, nso, 1.0, Cp[0], nmo, Vmizp[0], nocc, 0.0, Vpip[0], nocc);
            psio_->write(PSIF_HESS, "Vpi^A", (char*)Vpip[0], static_cast<size_t>(nmo) * nocc * sizeof(double), next_Vpi,
                         &next_Vpi);
        }
    }

    // => Jpi/Kpi <= //
    {
        // Figure out DFT functional info
        double Kscale = functional_->x_alpha();
        if (functional_->is_x_lrc()) throw PSIEXCEPTION("Hessians for LRC functionals are not implemented yet.");

        size_t memory = 0.9 * memory_ / 8L;
        size_t max_a = memory / (3L * nso * nso);
        max_a = (max_a > 3 * natom ? 3 * natom : max_a);

        int natom = basisset_->molecule()->natom();

        std::vector<SharedMatrix> dGmats;
        std::vector<double**> pdG(3 * natom);
        std::vector<bool> pert_incore(3 * natom);
        for (int a = 0; a < max_a; ++a)
            dGmats.push_back(std::make_shared<Matrix>("G derivative contribution", nso, nso));

        auto Gpi = std::make_shared<Matrix>("MO G Deriv", nmo, nocc);
        double** pGpi = Gpi->pointer();
        psio_address next_Gpi = PSIO_ZERO;
        // Write some junk for now, to set the sizing for PSIO
        for (int pert = 0; pert < 3 * natom; ++pert) {
            psio_->write(PSIF_HESS, "Gpi^A", (char*)pGpi[0], static_cast<size_t>(nmo) * nocc * sizeof(double), next_Gpi,
                         &next_Gpi);
        }
        next_Gpi = PSIO_ZERO;

        if (options_.get_str("SCF_TYPE").find("DF") != std::string::npos) {
            /*
             *  The DF algorithm
             */
            std::shared_ptr<BasisSet> auxiliary_ = get_basisset("DF_BASIS_SCF");

            auto Pmnfactory =
                std::make_shared<IntegralFactory>(auxiliary_, BasisSet::zero_ao_basis_set(), basisset_, basisset_);
            auto PQfactory = std::make_shared<IntegralFactory>(auxiliary_, BasisSet::zero_ao_basis_set(), auxiliary_,
                                                               BasisSet::zero_ao_basis_set());
            std::shared_ptr<TwoBodyAOInt> Pmnint(Pmnfactory->eri(2));
            std::shared_ptr<TwoBodyAOInt> PQint(PQfactory->eri(2));
            int np = auxiliary_->nbf();
            int nso = basisset_->nbf();
            int nauxshell = auxiliary_->nshell();
            int nshell = basisset_->nshell();
            int maxp = auxiliary_->max_function_per_shell();

            auto Amn = std::make_shared<Matrix>("(A|mn)", np, nso * nso);
            auto Ami = std::make_shared<Matrix>("(A|mi)", np, nso * nocc);
            auto Aij = std::make_shared<Matrix>("(A|ij)", np, nocc * nocc);
            auto Bmn = std::make_shared<Matrix>("Minv[B][A] (A|mn)", np, nso * nso);
            auto Tmn = std::make_shared<Matrix>("Tmn", np, nso * nso);
            auto TempP = std::make_shared<Matrix>("Temp[P]", 9, maxp);
            auto TempPmn = std::make_shared<Matrix>("Temp[P][mn]", maxp, nso * nso);
            auto c = std::make_shared<Vector>("c[A] = (mn|A) D[m][n]", np);
            auto d = std::make_shared<Vector>("d[A] = Minv[A][B] C[B]", np);
            double** Amnp = Amn->pointer();
            double** Amip = Ami->pointer();
            double** Aijp = Aij->pointer();
            double** Bmnp = Bmn->pointer();
            double** pTmn = Tmn->pointer();
            double** pTempP = TempP->pointer();
            double** pTmpPmn = TempPmn->pointer();
            double* cp = c->pointer();
            double* dp = d->pointer();

            // This probably shouldn't be recomputed here; we already needed it to get the
            // second derivative integrals.  One fine day, this should be refactored.
            auto metric = std::make_shared<FittingMetric>(auxiliary_, true);
            metric->form_full_eig_inverse(options_.get_double("DF_FITTING_CONDITION"));
            SharedMatrix PQ = metric->get_metric();
            double** PQp = PQ->pointer();

            // Same applies to these terms.  There are already hooks to compute c and d vectors, and store them on disk,
            // so we should make better use of those intermediates between the second derivative integrals and these
            // first derivative terms needed for the Fock matrix derivatives.
            for (int P = 0; P < nauxshell; ++P) {
                int nP = auxiliary_->shell(P).nfunction();
                int oP = auxiliary_->shell(P).function_index();
                for (int M = 0; M < nshell; ++M) {
                    int nM = basisset_->shell(M).nfunction();
                    int oM = basisset_->shell(M).function_index();
                    for (int N = 0; N < nshell; ++N) {
                        int nN = basisset_->shell(N).nfunction();
                        int oN = basisset_->shell(N).function_index();

                        Pmnint->compute_shell(P, 0, M, N);
                        const double* buffer = Pmnint->buffer();

                        for (int p = oP; p < oP+nP; p++) {
                            for (int m = oM; m < oM+nM; m++) {
                                for (int n = oN; n < oN+nN; n++) {
                                    Amnp[p][m*nso+n] = (*buffer++);
                                }
                            }
                        }
                    }
                }
            }
            // c[A] = (A|mn) D[m][n]
            C_DGEMV('N', np, nso*(size_t)nso, 1.0, Amnp[0], nso*(size_t)nso, Dap[0], 1, 0.0, cp, 1);
            // (A|mj) = (A|mn) C[n][j]
            C_DGEMM('N','N',np*(size_t)nso,nocc,nso,1.0,Amnp[0],nso,Cop[0],nocc,0.0,Amip[0],nocc);
            // (A|ij) = (A|mj) C[m][i]
            #pragma omp parallel for
            for (int p = 0; p < np; p++) {
                C_DGEMM('T','N',nocc,nocc,nso,1.0,Amip[p],nocc,Cop[0],nocc,0.0,&Aijp[0][p * (size_t) nocc * nocc],nocc);
            }

            // d[A] = Minv[A][B] c[B]  (factor of 2, to account for RHF)
            C_DGEMV('n', np, np, 2.0, PQp[0], np, cp, 1, 0.0, dp, 1);

            // B[B][m,n] = Minv[A][B] (A|mn)
            C_DGEMM('n', 'n', np, nso * nso, np, 1.0, PQp[0], np, Amnp[0], nso * nso, 0.0, Bmnp[0], nso * nso);

            // T[p][m,n] = B[p][r,n] D[m,r]
#pragma omp parallel for
            for (int p = 0; p < np; ++p)
                C_DGEMM('t', 'n', nso, nso, nso, 1.0, Dap[0], nso, Bmnp[p], nso, 0.0, pTmn[p], nso);

            for (int A = 0; A < 3 * natom; A += max_a) {
                int nA = (A + max_a >= 3 * natom ? 3 * natom - A : max_a);

                // Keep track of which centers are loaded into memory, so we know when to skip
                std::fill(pert_incore.begin(), pert_incore.end(), false);
                std::fill(pdG.begin(), pdG.end(), (double**)nullptr);
                for (int a = 0; a < nA; a++) {
                    pert_incore[floor((A + a) / 3.0)] = true;
                    pdG[A + a] = dGmats[a]->pointer();
                    dGmats[a]->zero();
                }

                for (int P = 0; P < nauxshell; ++P) {
                    int nP = auxiliary_->shell(P).nfunction();
                    int oP = auxiliary_->shell(P).function_index();
                    int Pcenter = auxiliary_->shell(P).ncenter();
                    int Pncart = auxiliary_->shell(P).ncartesian();
                    int Px = 3 * Pcenter + 0;
                    int Py = 3 * Pcenter + 1;
                    int Pz = 3 * Pcenter + 2;
                    for (int Q = 0; Q < nauxshell; ++Q) {
                        int nQ = auxiliary_->shell(Q).nfunction();
                        int oQ = auxiliary_->shell(Q).function_index();
                        int Qcenter = auxiliary_->shell(Q).ncenter();
                        int Qncart = auxiliary_->shell(Q).ncartesian();
                        int Qx = 3 * Qcenter + 0;
                        int Qy = 3 * Qcenter + 1;
                        int Qz = 3 * Qcenter + 2;


                        if (!pert_incore[Pcenter] && !pert_incore[Qcenter]) continue;

                        PQint->compute_shell_deriv1(P, 0, Q, 0);
                        const auto& buffers = PQint->buffers();
                        const double* PxBuf = buffers[0];
                        const double* PyBuf = buffers[1];
                        const double* PzBuf = buffers[2];
                        const double* QxBuf = buffers[3];
                        const double* QyBuf = buffers[4];
                        const double* QzBuf = buffers[5];

                        if (pert_incore[Pcenter]) {
                            // J terms
                            // Px
                            C_DGEMV('n', nP, nQ, 1.0, const_cast<double*>(PxBuf), nQ, &dp[oQ], 1, 0.0, pTempP[0], 1);
                            C_DGEMV('t', nP, nso * nso, -1.0, Bmnp[oP], nso * nso, pTempP[0], 1, 1.0, pdG[Px][0], 1);
                            // Py
                            C_DGEMV('n', nP, nQ, 1.0, const_cast<double*>(PyBuf), nQ, &dp[oQ], 1, 0.0, pTempP[0], 1);
                            C_DGEMV('t', nP, nso * nso, -1.0, Bmnp[oP], nso * nso, pTempP[0], 1, 1.0, pdG[Py][0], 1);
                            // Pz
                            C_DGEMV('n', nP, nQ, 1.0, const_cast<double*>(PzBuf), nQ, &dp[oQ], 1, 0.0, pTempP[0], 1);
                            C_DGEMV('t', nP, nso * nso, -1.0, Bmnp[oP], nso * nso, pTempP[0], 1, 1.0, pdG[Pz][0], 1);

                            if (Kscale) {
                                // K terms
                                // Px
                                C_DGEMM('n', 'n', nP, nso * nso, nQ, 1.0, const_cast<double*>(PxBuf), nQ, pTmn[oQ],
                                        nso * nso, 0.0, pTmpPmn[0], nso * nso);
                                for (int p = 0; p < nP; ++p)
                                    C_DGEMM('N', 'N', nso, nso, nso, Kscale, Bmnp[p + oP], nso, pTmpPmn[p], nso, 1.0,
                                            pdG[Px][0], nso);
                                // Py
                                C_DGEMM('n', 'n', nP, nso * nso, nQ, 1.0, const_cast<double*>(PyBuf), nQ, pTmn[oQ],
                                        nso * nso, 0.0, pTmpPmn[0], nso * nso);
                                for (int p = 0; p < nP; ++p)
                                    C_DGEMM('N', 'N', nso, nso, nso, Kscale, Bmnp[p + oP], nso, pTmpPmn[p], nso, 1.0,
                                            pdG[Py][0], nso);
                                // Pz
                                C_DGEMM('n', 'n', nP, nso * nso, nQ, 1.0, const_cast<double*>(PzBuf), nQ, pTmn[oQ],
                                        nso * nso, 0.0, pTmpPmn[0], nso * nso);
                                for (int p = 0; p < nP; ++p)
                                    C_DGEMM('N', 'N', nso, nso, nso, Kscale, Bmnp[p + oP], nso, pTmpPmn[p], nso, 1.0,
                                            pdG[Pz][0], nso);
                            }
                        }
                        if (pert_incore[Qcenter]) {
                            // J terms
                            // Qx
                            C_DGEMV('n', nP, nQ, 1.0, const_cast<double*>(QxBuf), nQ, &dp[oQ], 1, 0.0, pTempP[0], 1);
                            C_DGEMV('t', nP, nso * nso, -1.0, Bmnp[oP], nso * nso, pTempP[0], 1, 1.0, pdG[Qx][0], 1);
                            // Qy
                            C_DGEMV('n', nP, nQ, 1.0, const_cast<double*>(QyBuf), nQ, &dp[oQ], 1, 0.0, pTempP[0], 1);
                            C_DGEMV('t', nP, nso * nso, -1.0, Bmnp[oP], nso * nso, pTempP[0], 1, 1.0, pdG[Qy][0], 1);
                            // Qz
                            C_DGEMV('n', nP, nQ, 1.0, const_cast<double*>(QzBuf), nQ, &dp[oQ], 1, 0.0, pTempP[0], 1);
                            C_DGEMV('t', nP, nso * nso, -1.0, Bmnp[oP], nso * nso, pTempP[0], 1, 1.0, pdG[Qz][0], 1);

                            if (Kscale) {
                                // K terms
                                // Qx
                                C_DGEMM('n', 'n', nP, nso * nso, nQ, 1.0, const_cast<double*>(QxBuf), nQ, pTmn[oQ],
                                        nso * nso, 0.0, pTmpPmn[0], nso * nso);
                                for (int p = 0; p < nP; ++p)
                                    C_DGEMM('N', 'N', nso, nso, nso, Kscale, Bmnp[p + oP], nso, pTmpPmn[p], nso, 1.0,
                                            pdG[Qx][0], nso);
                                // Qy
                                C_DGEMM('n', 'n', nP, nso * nso, nQ, 1.0, const_cast<double*>(QyBuf), nQ, pTmn[oQ],
                                        nso * nso, 0.0, pTmpPmn[0], nso * nso);
                                for (int p = 0; p < nP; ++p)
                                    C_DGEMM('N', 'N', nso, nso, nso, Kscale, Bmnp[p + oP], nso, pTmpPmn[p], nso, 1.0,
                                            pdG[Qy][0], nso);
                                // Qz
                                C_DGEMM('n', 'n', nP, nso * nso, nQ, 1.0, const_cast<double*>(QzBuf), nQ, pTmn[oQ],
                                        nso * nso, 0.0, pTmpPmn[0], nso * nso);
                                for (int p = 0; p < nP; ++p)
                                    C_DGEMM('N', 'N', nso, nso, nso, Kscale, Bmnp[p + oP], nso, pTmpPmn[p], nso, 1.0,
                                            pdG[Qz][0], nso);
                            }
                        }
                    }
                }

                for (int P = 0; P < nauxshell; ++P) {
                    int nP = auxiliary_->shell(P).nfunction();
                    int oP = auxiliary_->shell(P).function_index();
                    int Pcenter = auxiliary_->shell(P).ncenter();
                    int Pncart = auxiliary_->shell(P).ncartesian();
                    int Px = 3 * Pcenter + 0;
                    int Py = 3 * Pcenter + 1;
                    int Pz = 3 * Pcenter + 2;
                    for (int M = 0; M < nshell; ++M) {
                        int nM = basisset_->shell(M).nfunction();
                        int oM = basisset_->shell(M).function_index();
                        int Mcenter = basisset_->shell(M).ncenter();
                        int Mncart = basisset_->shell(M).ncartesian();
                        int mx = 3 * Mcenter + 0;
                        int my = 3 * Mcenter + 1;
                        int mz = 3 * Mcenter + 2;
                        for (int N = 0; N < nshell; ++N) {
                            int nN = basisset_->shell(N).nfunction();
                            int oN = basisset_->shell(N).function_index();
                            int Ncenter = basisset_->shell(N).ncenter();
                            int Nncart = basisset_->shell(N).ncartesian();
                            int nx = 3 * Ncenter + 0;
                            int ny = 3 * Ncenter + 1;
                            int nz = 3 * Ncenter + 2;


                            if (!pert_incore[Pcenter] && !pert_incore[Mcenter] && !pert_incore[Ncenter]) continue;

                            Pmnint->compute_shell_deriv1(P, 0, M, N);
                            const auto& buffers = Pmnint->buffers();
                            const double* PxBuf = buffers[0];
                            const double* PyBuf = buffers[1];
                            const double* PzBuf = buffers[2];
                            const double* mxBuf = buffers[3];
                            const double* myBuf = buffers[4];
                            const double* mzBuf = buffers[5];
                            const double* nxBuf = buffers[6];
                            const double* nyBuf = buffers[7];
                            const double* nzBuf = buffers[8];

                            /*
                             * J terms have 2 contributions:
                             *      F^x[m][n] <- (P|mn)^x d[P]
                             * and
                             *      F^x[r][s] <- D[m][n] (P|mn)^x B[P][r,s]
                             * The second term factorizes into...
                             * ... Temp[P] = D[m][n] (P|mn)^x ...
                             * ... and then F^x[r][s] <- Temp[P] B[P][r,s]  (factor of 2 for RHF)
                             */

                            size_t delta = 0L;
                            for (int p = 0; p < nP; ++p) {
                                double valPx = 0, valPy = 0, valPz = 0;
                                double valmx = 0, valmy = 0, valmz = 0;
                                double valnx = 0, valny = 0, valnz = 0;
                                for (int m = oM; m < nM + oM; ++m) {
                                    for (int n = oN; n < nN + oN; ++n) {
                                        valPx += Dap[m][n] * PxBuf[delta];
                                        valPy += Dap[m][n] * PyBuf[delta];
                                        valPz += Dap[m][n] * PzBuf[delta];
                                        valmx += Dap[m][n] * mxBuf[delta];
                                        valmy += Dap[m][n] * myBuf[delta];
                                        valmz += Dap[m][n] * mzBuf[delta];
                                        valnx += Dap[m][n] * nxBuf[delta];
                                        valny += Dap[m][n] * nyBuf[delta];
                                        valnz += Dap[m][n] * nzBuf[delta];
                                        ++delta;
                                    }
                                }
                                pTempP[0][p] = valPx;
                                pTempP[1][p] = valPy;
                                pTempP[2][p] = valPz;
                                pTempP[3][p] = valmx;
                                pTempP[4][p] = valmy;
                                pTempP[5][p] = valmz;
                                pTempP[6][p] = valnx;
                                pTempP[7][p] = valny;
                                pTempP[8][p] = valnz;
                            }

                            if (pert_incore[Pcenter]) {
                                // J Terms
                                size_t delta = 0L;
                                for (int p = oP; p < oP + nP; ++p) {
                                    for (int m = oM; m < nM + oM; ++m) {
                                        for (int n = oN; n < nN + oN; ++n) {
                                            pdG[Px][m][n] += PxBuf[delta] * dp[p];
                                            pdG[Py][m][n] += PyBuf[delta] * dp[p];
                                            pdG[Pz][m][n] += PzBuf[delta] * dp[p];
                                            ++delta;
                                        }
                                    }
                                }
                                C_DGEMV('t', nP, nso * nso, 2.0, Bmnp[oP], nso * nso, pTempP[0], 1, 1.0, pdG[Px][0], 1);
                                C_DGEMV('t', nP, nso * nso, 2.0, Bmnp[oP], nso * nso, pTempP[1], 1, 1.0, pdG[Py][0], 1);
                                C_DGEMV('t', nP, nso * nso, 2.0, Bmnp[oP], nso * nso, pTempP[2], 1, 1.0, pdG[Pz][0], 1);
                                // K Terms
                                if (Kscale) {
                                    for (int p = 0; p < nP; ++p)
                                        C_DGEMM('T', 'N', nN, nso, nM, -2 * Kscale,
                                                const_cast<double*>(PxBuf) + p * nM * nN, nN, &pTmn[oP + p][oM * nso],
                                                nso, 1.0, pdG[Px][oN], nso);
                                    for (int p = 0; p < nP; ++p)
                                        C_DGEMM('T', 'N', nN, nso, nM, -2 * Kscale,
                                                const_cast<double*>(PyBuf) + p * nM * nN, nN, &pTmn[oP + p][oM * nso],
                                                nso, 1.0, pdG[Py][oN], nso);
                                    for (int p = 0; p < nP; ++p)
                                        C_DGEMM('T', 'N', nN, nso, nM, -2 * Kscale,
                                                const_cast<double*>(PzBuf) + p * nM * nN, nN, &pTmn[oP + p][oM * nso],
                                                nso, 1.0, pdG[Pz][oN], nso);
                                }
                            }
                            if (pert_incore[Mcenter]) {
                                // J Terms
                                size_t delta = 0L;
                                for (int p = oP; p < oP + nP; ++p) {
                                    for (int m = oM; m < nM + oM; ++m) {
                                        for (int n = oN; n < nN + oN; ++n) {
                                            pdG[mx][m][n] += mxBuf[delta] * dp[p];
                                            pdG[my][m][n] += myBuf[delta] * dp[p];
                                            pdG[mz][m][n] += mzBuf[delta] * dp[p];
                                            ++delta;
                                        }
                                    }
                                }
                                C_DGEMV('t', nP, nso * nso, 2.0, Bmnp[oP], nso * nso, pTempP[3], 1, 1.0, pdG[mx][0], 1);
                                C_DGEMV('t', nP, nso * nso, 2.0, Bmnp[oP], nso * nso, pTempP[4], 1, 1.0, pdG[my][0], 1);
                                C_DGEMV('t', nP, nso * nso, 2.0, Bmnp[oP], nso * nso, pTempP[5], 1, 1.0, pdG[mz][0], 1);
                                // K Terms
                                if (Kscale) {
                                    for (int p = 0; p < nP; ++p)
                                        C_DGEMM('T', 'N', nN, nso, nM, -2 * Kscale,
                                                const_cast<double*>(mxBuf) + p * nM * nN, nN, &pTmn[oP + p][oM * nso],
                                                nso, 1.0, pdG[mx][oN], nso);
                                    for (int p = 0; p < nP; ++p)
                                        C_DGEMM('T', 'N', nN, nso, nM, -2 * Kscale,
                                                const_cast<double*>(myBuf) + p * nM * nN, nN, &pTmn[oP + p][oM * nso],
                                                nso, 1.0, pdG[my][oN], nso);
                                    for (int p = 0; p < nP; ++p)
                                        C_DGEMM('T', 'N', nN, nso, nM, -2 * Kscale,
                                                const_cast<double*>(mzBuf) + p * nM * nN, nN, &pTmn[oP + p][oM * nso],
                                                nso, 1.0, pdG[mz][oN], nso);
                                }
                            }
                            if (pert_incore[Ncenter]) {
                                // J Terms
                                size_t delta = 0L;
                                for (int p = oP; p < oP + nP; ++p) {
                                    for (int m = oM; m < nM + oM; ++m) {
                                        for (int n = oN; n < nN + oN; ++n) {
                                            pdG[nx][m][n] += nxBuf[delta] * dp[p];
                                            pdG[ny][m][n] += nyBuf[delta] * dp[p];
                                            pdG[nz][m][n] += nzBuf[delta] * dp[p];
                                            ++delta;
                                        }
                                    }
                                }
                                C_DGEMV('t', nP, nso * nso, 2.0, Bmnp[oP], nso * nso, pTempP[6], 1, 1.0, pdG[nx][0], 1);
                                C_DGEMV('t', nP, nso * nso, 2.0, Bmnp[oP], nso * nso, pTempP[7], 1, 1.0, pdG[ny][0], 1);
                                C_DGEMV('t', nP, nso * nso, 2.0, Bmnp[oP], nso * nso, pTempP[8], 1, 1.0, pdG[nz][0], 1);
                                // K Terms
                                if (Kscale) {
                                    for (int p = 0; p < nP; ++p)
                                        C_DGEMM('T', 'N', nN, nso, nM, -2 * Kscale,
                                                const_cast<double*>(nxBuf) + p * nM * nN, nN, &pTmn[oP + p][oM * nso],
                                                nso, 1.0, pdG[nx][oN], nso);
                                    for (int p = 0; p < nP; ++p)
                                        C_DGEMM('T', 'N', nN, nso, nM, -2 * Kscale,
                                                const_cast<double*>(nyBuf) + p * nM * nN, nN, &pTmn[oP + p][oM * nso],
                                                nso, 1.0, pdG[ny][oN], nso);
                                    for (int p = 0; p < nP; ++p)
                                        C_DGEMM('T', 'N', nN, nso, nM, -2 * Kscale,
                                                const_cast<double*>(nzBuf) + p * nM * nN, nN, &pTmn[oP + p][oM * nso],
                                                nso, 1.0, pdG[nz][oN], nso);
                                }
                            }
                        }
                    }
                }

                for (int a = 0; a < nA; ++a) {
                    // Symmetrize the derivative Fock contributions
                    SharedMatrix G = dGmats[a];
                    G->add(G->transpose());
                    Gpi->transform(C, G, Cocc);
                    Gpi->scale(0.5);
                    psio_->write(PSIF_HESS, "Gpi^A", (char*)pGpi[0], static_cast<size_t>(nmo) * nocc * sizeof(double),
                                 next_Gpi, &next_Gpi);
                }

            }  // End loop over A batches

        } else {
            /*
             * The conventional integral algorithm
             */
#ifdef USING_BrianQC
            if (brianEnable) {
                brianBool computeCoulomb = BRIAN_TRUE;
                brianBool computeExchange = BRIAN_TRUE;
                
                brianInt maxSegmentSize;
                brianCPHFMaxSegmentSize(&brianCookie, &maxSegmentSize);
                brianInt maxSegmentAtomCount = maxSegmentSize / 3;
                brianInt segmentAtomCount = -1;
                brianInt segmentAtomIndexStart = -1;
                
                std::vector<std::shared_ptr<Matrix>> Jmnx(maxSegmentAtomCount);
                std::vector<std::shared_ptr<Matrix>> Jmny(maxSegmentAtomCount);
                std::vector<std::shared_ptr<Matrix>> Jmnz(maxSegmentAtomCount);
                std::vector<std::shared_ptr<Matrix>> Kmnx(maxSegmentAtomCount);
                std::vector<std::shared_ptr<Matrix>> Kmny(maxSegmentAtomCount);
                std::vector<std::shared_ptr<Matrix>> Kmnz(maxSegmentAtomCount);
                std::vector<double*> Jmn(maxSegmentAtomCount * 3);
                std::vector<double*> Kmn(maxSegmentAtomCount * 3);
                for (int i = 0; i < maxSegmentAtomCount; i++) {
                    Jmnx[i] = std::make_shared<Matrix>("Jmnx", nso, nso);
                    Jmny[i] = std::make_shared<Matrix>("Jmny", nso, nso);
                    Jmnz[i] = std::make_shared<Matrix>("Jmnz", nso, nso);
                    Kmnx[i] = std::make_shared<Matrix>("Kmnx", nso, nso);
                    Kmny[i] = std::make_shared<Matrix>("Kmny", nso, nso);
                    Kmnz[i] = std::make_shared<Matrix>("Kmnz", nso, nso);
                    Jmn[i * 3 + 0] = Jmnx[i]->get_pointer();
                    Jmn[i * 3 + 1] = Jmny[i]->get_pointer();
                    Jmn[i * 3 + 2] = Jmnz[i]->get_pointer();
                    Kmn[i * 3 + 0] = Kmnx[i]->get_pointer();
                    Kmn[i * 3 + 1] = Kmny[i]->get_pointer();
                    Kmn[i * 3 + 2] = Kmnz[i]->get_pointer();
                }
                
                for (int A = 0; A < natom; A++) {
                    if (segmentAtomCount < 0 || A < segmentAtomIndexStart || A >= (segmentAtomIndexStart + segmentAtomCount)) {
                        segmentAtomIndexStart = A;
                        segmentAtomCount = (segmentAtomIndexStart + maxSegmentAtomCount > natom) ? (natom - segmentAtomIndexStart) : maxSegmentAtomCount;
                        
                        brianInt integralType = BRIAN_INTEGRAL_TYPE_NUCLEAR;
                        brianCPHFBuildRepulsionDeriv(&brianCookie, &computeCoulomb, &computeExchange, &segmentAtomCount, &segmentAtomIndexStart, Dt->get_pointer(), nullptr, Jmn.data(), Kmn.data(), nullptr);
                    }
                    
                    Jmnx[A - segmentAtomIndexStart]->subtract(Kmnx[A - segmentAtomIndexStart]);
                    Jmny[A - segmentAtomIndexStart]->subtract(Kmny[A - segmentAtomIndexStart]);
                    Jmnz[A - segmentAtomIndexStart]->subtract(Kmnz[A - segmentAtomIndexStart]);
                    
                    Gpi->transform(C, Jmnx[A - segmentAtomIndexStart], Cocc);
                    psio_->write(PSIF_HESS,"Gpi^A",(char*)pGpi[0],nmo * nocc * sizeof(double),next_Gpi,&next_Gpi);
                    Gpi->transform(C, Jmny[A - segmentAtomIndexStart], Cocc);
                    psio_->write(PSIF_HESS,"Gpi^A",(char*)pGpi[0],nmo * nocc * sizeof(double),next_Gpi,&next_Gpi);
                    Gpi->transform(C, Jmnz[A - segmentAtomIndexStart], Cocc);
                    psio_->write(PSIF_HESS,"Gpi^A",(char*)pGpi[0],nmo * nocc * sizeof(double),next_Gpi,&next_Gpi);
                }
            } else {
#endif
            std::shared_ptr<TwoBodyAOInt> ints(integral_->eri(1));

            const std::vector<std::pair<int, int>>& shell_pairs = ints->shell_pairs();
            size_t npairs = shell_pairs.size();
            size_t npairs2 = npairs * npairs;

            std::vector<std::mutex> mutexes(3 * natom);
            for (int A = 0; A < 3 * natom; A += max_a) {
                int nA = (A + max_a >= 3 * natom ? 3 * natom - A : max_a);

                // Keep track of which centers are loaded into memory, so we know when to skip
                std::fill(pert_incore.begin(), pert_incore.end(), false);
                std::fill(pdG.begin(), pdG.end(), (double**)nullptr);
                for (int a = 0; a < nA; a++) {
                    pert_incore[floor((A + a) / 3.0)] = true;
                    pdG[A + a] = dGmats[a]->pointer();
                    dGmats[a]->zero();
                }

                int maxfunc = basisset_->max_function_per_shell();
                int maxnpair = maxfunc * maxfunc;
                std::vector<std::shared_ptr<TwoBodyAOInt>> ints;
                // clang-format off
                    // Ordering of per-thread scratch buffers:
                    //    Jpq     |    Jrs    |    Kpr    |    Kqs    |    Kps    |    Kqr
                    //  ---------------------------------------------------------------------
                    //  x   y   z | x   y   z | x   y   z | x   y   z | x   y   z | x   y   z
                    //  0   1   2 | 3   4   5 | 6   7   8 | 9  10  11 | 12 13  14 |15  16  17
                // clang-format on
                int nbuffers = (Kscale != 0.0 ? 6 : 2) * 3;
                int nthreads = Process::environment.get_n_threads();
                std::vector<std::vector<std::vector<double>>> thread_temps;
                for (int thread = 0; thread < nthreads; thread++) {
                    ints.push_back(std::shared_ptr<TwoBodyAOInt>(integral_->eri(1)));
                    std::vector<std::vector<double>> temps(4);
                    for (auto& t : temps) t.resize((size_t)nbuffers * maxnpair);
                    thread_temps.push_back(temps);
                }
                size_t computed_shells = 0L;
                // shell pair blocks
                auto blocksPQ = ints[0]->get_blocks12();
                auto blocksRS = ints[0]->get_blocks34();
                bool use_batching = blocksPQ != blocksRS;
#pragma omp parallel for schedule(dynamic) num_threads(nthreads)
                // loop over all the blocks of (P>=Q|
                for (size_t blockPQ_idx = 0; blockPQ_idx < blocksPQ.size(); blockPQ_idx++) {
                    const auto& blockPQ = blocksPQ[blockPQ_idx];
#ifdef _OPENMP
                    const int rank = omp_get_thread_num();
#else
                    const int rank = 0;
#endif
                    auto& buffers = ints[rank]->buffers();
                    // loop over all the blocks of |R>=S)
                    int loop_start = use_batching ? 0 : blockPQ_idx;
                    for (int blockRS_idx = loop_start; blockRS_idx < blocksRS.size(); ++blockRS_idx) {
                        const auto& blockRS = blocksRS[blockRS_idx];

                        // Check to see if any of these integrals could contribute
                        bool contributes = false;
                        for (const auto& pairPQ : blockPQ) {
                            const auto& Pcenter = basisset_->shell(pairPQ.first).ncenter();
                            const auto& Qcenter = basisset_->shell(pairPQ.second).ncenter();
                            if (pert_incore[Pcenter] || pert_incore[Qcenter]) {
                                contributes = true;
                                break;
                            }
                        }
                        if (!contributes) {
                            for (const auto& pairRS : blockRS) {
                                const auto& Rcenter = basisset_->shell(pairRS.first).ncenter();
                                const auto& Scenter = basisset_->shell(pairRS.second).ncenter();
                                if (pert_incore[Rcenter] || pert_incore[Scenter]) {
                                    contributes = true;
                                    break;
                                }
                            }
                        }
                        if (!contributes) continue;

                        // This is where we want to screen with density and schwarz-like screening
                        // compute the integrals and continue if none were computed
                        ints[rank]->compute_shell_blocks_deriv1(blockPQ_idx, blockRS_idx);

                        const double* pAx = buffers[0];
                        const double* pAy = buffers[1];
                        const double* pAz = buffers[2];
                        const double* pBx = buffers[3];
                        const double* pBy = buffers[4];
                        const double* pBz = buffers[5];
                        const double* pCx = buffers[6];
                        const double* pCy = buffers[7];
                        const double* pCz = buffers[8];
                        const double* pDx = buffers[9];
                        const double* pDy = buffers[10];
                        const double* pDz = buffers[11];
                        std::array<double*, 4> pTemps;
                        for (int i = 0; i < 4; ++i) pTemps[i] = thread_temps[rank][i].data();
                        // Loop over all of the P,Q,R,S shells within the blocks.  We have P>=Q, R>=S and PQ<=RS.
                        for (const auto& pairPQ : blockPQ) {
                            const auto& P = pairPQ.first;
                            const auto& Q = pairPQ.second;
                            const auto& Pshell = basisset_->shell(P);
                            const auto& Qshell = basisset_->shell(Q);
                            const auto& Psize = Pshell.nfunction();
                            const auto& Qsize = Qshell.nfunction();
                            const auto& Poff = Pshell.function_index();
                            const auto& Qoff = Qshell.function_index();
                            const auto& Pcenter = Pshell.ncenter();
                            const auto& Qcenter = Qshell.ncenter();

                            for (const auto& pairRS : blockRS) {
                                const auto& R = pairRS.first;
                                const auto& S = pairRS.second;
                                const auto& Rshell = basisset_->shell(R);
                                const auto& Sshell = basisset_->shell(S);
                                const auto& Rsize = Rshell.nfunction();
                                const auto& Ssize = Sshell.nfunction();
                                const auto& Roff = Rshell.function_index();
                                const auto& Soff = Sshell.function_index();
                                const auto& Rcenter = Rshell.ncenter();
                                const auto& Scenter = Sshell.ncenter();

                                size_t block_size = (size_t)Psize * Qsize * Rsize * Ssize;
                                // When there are chunks of shellpairs in RS, we need to make sure
                                // we filter out redundant combinations.  This should probably be done
                                // by having a block of RS generated for each PQ at list build time.
                                if (use_batching && ((P > R) || (P == R && Q > S))) {
                                    pAx += block_size;
                                    pAy += block_size;
                                    pAz += block_size;
                                    pBx += block_size;
                                    pBy += block_size;
                                    pBz += block_size;
                                    pCx += block_size;
                                    pCy += block_size;
                                    pCz += block_size;
                                    pDx += block_size;
                                    pDy += block_size;
                                    pDz += block_size;
                                    continue;
                                }

                                double prefactor = 8.0;
                                if (P == Q) prefactor *= 0.5;
                                if (R == S) prefactor *= 0.5;
                                if (P == R && Q == S) prefactor *= 0.5;

                                double Dpq, Drs, Dpr, Dqs, Dps, Dqr;
                                size_t delta;
                                double Ax, Ay, Az;
                                double Bx, By, Bz;
                                double Cx, Cy, Cz;
                                double Dx, Dy, Dz;

                                // => Coulomb Term <= //
                                if (pert_incore[Pcenter]) {
                                    std::fill_n(pTemps[0] + 0 * maxnpair, Psize * Qsize, 0.0);
                                    std::fill_n(pTemps[0] + 1 * maxnpair, Psize * Qsize, 0.0);
                                    std::fill_n(pTemps[0] + 2 * maxnpair, Psize * Qsize, 0.0);
                                    std::fill_n(pTemps[0] + 3 * maxnpair, Rsize * Ssize, 0.0);
                                    std::fill_n(pTemps[0] + 4 * maxnpair, Rsize * Ssize, 0.0);
                                    std::fill_n(pTemps[0] + 5 * maxnpair, Rsize * Ssize, 0.0);
                                }
                                if (pert_incore[Qcenter]) {
                                    std::fill_n(pTemps[1] + 0 * maxnpair, Psize * Qsize, 0.0);
                                    std::fill_n(pTemps[1] + 1 * maxnpair, Psize * Qsize, 0.0);
                                    std::fill_n(pTemps[1] + 2 * maxnpair, Psize * Qsize, 0.0);
                                    std::fill_n(pTemps[1] + 3 * maxnpair, Rsize * Ssize, 0.0);
                                    std::fill_n(pTemps[1] + 4 * maxnpair, Rsize * Ssize, 0.0);
                                    std::fill_n(pTemps[1] + 5 * maxnpair, Rsize * Ssize, 0.0);
                                }
                                if (pert_incore[Rcenter]) {
                                    std::fill_n(pTemps[2] + 0 * maxnpair, Psize * Qsize, 0.0);
                                    std::fill_n(pTemps[2] + 1 * maxnpair, Psize * Qsize, 0.0);
                                    std::fill_n(pTemps[2] + 2 * maxnpair, Psize * Qsize, 0.0);
                                    std::fill_n(pTemps[2] + 3 * maxnpair, Rsize * Ssize, 0.0);
                                    std::fill_n(pTemps[2] + 4 * maxnpair, Rsize * Ssize, 0.0);
                                    std::fill_n(pTemps[2] + 5 * maxnpair, Rsize * Ssize, 0.0);
                                }
                                if (pert_incore[Scenter]) {
                                    std::fill_n(pTemps[3] + 0 * maxnpair, Psize * Qsize, 0.0);
                                    std::fill_n(pTemps[3] + 1 * maxnpair, Psize * Qsize, 0.0);
                                    std::fill_n(pTemps[3] + 2 * maxnpair, Psize * Qsize, 0.0);
                                    std::fill_n(pTemps[3] + 3 * maxnpair, Rsize * Ssize, 0.0);
                                    std::fill_n(pTemps[3] + 4 * maxnpair, Rsize * Ssize, 0.0);
                                    std::fill_n(pTemps[3] + 5 * maxnpair, Rsize * Ssize, 0.0);
                                }
                                delta = 0L;
                                for (int p = 0; p < Psize; p++) {
                                    for (int q = 0; q < Qsize; q++) {
                                        for (int r = 0; r < Rsize; r++) {
                                            for (int s = 0; s < Ssize; s++) {
                                                Dpq = Dap[p + Poff][q + Qoff];
                                                Drs = Dap[r + Roff][s + Soff];
                                                Ax = prefactor * pAx[delta];
                                                Ay = prefactor * pAy[delta];
                                                Az = prefactor * pAz[delta];
                                                Bx = prefactor * pBx[delta];
                                                By = prefactor * pBy[delta];
                                                Bz = prefactor * pBz[delta];
                                                Cx = prefactor * pCx[delta];
                                                Cy = prefactor * pCy[delta];
                                                Cz = prefactor * pCz[delta];
                                                Dx = prefactor * pDx[delta];
                                                Dy = prefactor * pDy[delta];
                                                Dz = prefactor * pDz[delta];

                                                if (pert_incore[Pcenter]) {
                                                    pTemps[0][0 * maxnpair + p * Qsize + q] += Ax * Drs;
                                                    pTemps[0][1 * maxnpair + p * Qsize + q] += Ay * Drs;
                                                    pTemps[0][2 * maxnpair + p * Qsize + q] += Az * Drs;
                                                    pTemps[0][3 * maxnpair + r * Ssize + s] += Ax * Dpq;
                                                    pTemps[0][4 * maxnpair + r * Ssize + s] += Ay * Dpq;
                                                    pTemps[0][5 * maxnpair + r * Ssize + s] += Az * Dpq;
                                                }
                                                if (pert_incore[Qcenter]) {
                                                    pTemps[1][0 * maxnpair + p * Qsize + q] += Bx * Drs;
                                                    pTemps[1][1 * maxnpair + p * Qsize + q] += By * Drs;
                                                    pTemps[1][2 * maxnpair + p * Qsize + q] += Bz * Drs;
                                                    pTemps[1][3 * maxnpair + r * Ssize + s] += Bx * Dpq;
                                                    pTemps[1][4 * maxnpair + r * Ssize + s] += By * Dpq;
                                                    pTemps[1][5 * maxnpair + r * Ssize + s] += Bz * Dpq;
                                                }
                                                if (pert_incore[Rcenter]) {
                                                    pTemps[2][0 * maxnpair + p * Qsize + q] += Cx * Drs;
                                                    pTemps[2][1 * maxnpair + p * Qsize + q] += Cy * Drs;
                                                    pTemps[2][2 * maxnpair + p * Qsize + q] += Cz * Drs;
                                                    pTemps[2][3 * maxnpair + r * Ssize + s] += Cx * Dpq;
                                                    pTemps[2][4 * maxnpair + r * Ssize + s] += Cy * Dpq;
                                                    pTemps[2][5 * maxnpair + r * Ssize + s] += Cz * Dpq;
                                                }
                                                if (pert_incore[Scenter]) {
                                                    pTemps[3][0 * maxnpair + p * Qsize + q] += Dx * Drs;
                                                    pTemps[3][1 * maxnpair + p * Qsize + q] += Dy * Drs;
                                                    pTemps[3][2 * maxnpair + p * Qsize + q] += Dz * Drs;
                                                    pTemps[3][3 * maxnpair + r * Ssize + s] += Dx * Dpq;
                                                    pTemps[3][4 * maxnpair + r * Ssize + s] += Dy * Dpq;
                                                    pTemps[3][5 * maxnpair + r * Ssize + s] += Dz * Dpq;
                                                }
                                                delta++;
                                            }
                                        }
                                    }
                                }
                                // => Exchange Term <= //
                                if (Kscale) {
                                    if (pert_incore[Pcenter]) {
                                        std::fill_n(pTemps[0] + 6 * maxnpair, Psize * Rsize, 0.0);
                                        std::fill_n(pTemps[0] + 7 * maxnpair, Psize * Rsize, 0.0);
                                        std::fill_n(pTemps[0] + 8 * maxnpair, Psize * Rsize, 0.0);
                                        std::fill_n(pTemps[0] + 9 * maxnpair, Qsize * Ssize, 0.0);
                                        std::fill_n(pTemps[0] + 10 * maxnpair, Qsize * Ssize, 0.0);
                                        std::fill_n(pTemps[0] + 11 * maxnpair, Qsize * Ssize, 0.0);
                                        std::fill_n(pTemps[0] + 12 * maxnpair, Psize * Ssize, 0.0);
                                        std::fill_n(pTemps[0] + 13 * maxnpair, Psize * Ssize, 0.0);
                                        std::fill_n(pTemps[0] + 14 * maxnpair, Psize * Ssize, 0.0);
                                        std::fill_n(pTemps[0] + 15 * maxnpair, Qsize * Rsize, 0.0);
                                        std::fill_n(pTemps[0] + 16 * maxnpair, Qsize * Rsize, 0.0);
                                        std::fill_n(pTemps[0] + 17 * maxnpair, Qsize * Rsize, 0.0);
                                    }
                                    if (pert_incore[Qcenter]) {
                                        std::fill_n(pTemps[1] + 6 * maxnpair, Psize * Rsize, 0.0);
                                        std::fill_n(pTemps[1] + 7 * maxnpair, Psize * Rsize, 0.0);
                                        std::fill_n(pTemps[1] + 8 * maxnpair, Psize * Rsize, 0.0);
                                        std::fill_n(pTemps[1] + 9 * maxnpair, Qsize * Ssize, 0.0);
                                        std::fill_n(pTemps[1] + 10 * maxnpair, Qsize * Ssize, 0.0);
                                        std::fill_n(pTemps[1] + 11 * maxnpair, Qsize * Ssize, 0.0);
                                        std::fill_n(pTemps[1] + 12 * maxnpair, Psize * Ssize, 0.0);
                                        std::fill_n(pTemps[1] + 13 * maxnpair, Psize * Ssize, 0.0);
                                        std::fill_n(pTemps[1] + 14 * maxnpair, Psize * Ssize, 0.0);
                                        std::fill_n(pTemps[1] + 15 * maxnpair, Qsize * Rsize, 0.0);
                                        std::fill_n(pTemps[1] + 16 * maxnpair, Qsize * Rsize, 0.0);
                                        std::fill_n(pTemps[1] + 17 * maxnpair, Qsize * Rsize, 0.0);
                                    }
                                    if (pert_incore[Rcenter]) {
                                        std::fill_n(pTemps[2] + 6 * maxnpair, Psize * Rsize, 0.0);
                                        std::fill_n(pTemps[2] + 7 * maxnpair, Psize * Rsize, 0.0);
                                        std::fill_n(pTemps[2] + 8 * maxnpair, Psize * Rsize, 0.0);
                                        std::fill_n(pTemps[2] + 9 * maxnpair, Qsize * Ssize, 0.0);
                                        std::fill_n(pTemps[2] + 10 * maxnpair, Qsize * Ssize, 0.0);
                                        std::fill_n(pTemps[2] + 11 * maxnpair, Qsize * Ssize, 0.0);
                                        std::fill_n(pTemps[2] + 12 * maxnpair, Psize * Ssize, 0.0);
                                        std::fill_n(pTemps[2] + 13 * maxnpair, Psize * Ssize, 0.0);
                                        std::fill_n(pTemps[2] + 14 * maxnpair, Psize * Ssize, 0.0);
                                        std::fill_n(pTemps[2] + 15 * maxnpair, Qsize * Rsize, 0.0);
                                        std::fill_n(pTemps[2] + 16 * maxnpair, Qsize * Rsize, 0.0);
                                        std::fill_n(pTemps[2] + 17 * maxnpair, Qsize * Rsize, 0.0);
                                    }
                                    if (pert_incore[Scenter]) {
                                        std::fill_n(pTemps[3] + 6 * maxnpair, Psize * Rsize, 0.0);
                                        std::fill_n(pTemps[3] + 7 * maxnpair, Psize * Rsize, 0.0);
                                        std::fill_n(pTemps[3] + 8 * maxnpair, Psize * Rsize, 0.0);
                                        std::fill_n(pTemps[3] + 9 * maxnpair, Qsize * Ssize, 0.0);
                                        std::fill_n(pTemps[3] + 10 * maxnpair, Qsize * Ssize, 0.0);
                                        std::fill_n(pTemps[3] + 11 * maxnpair, Qsize * Ssize, 0.0);
                                        std::fill_n(pTemps[3] + 12 * maxnpair, Psize * Ssize, 0.0);
                                        std::fill_n(pTemps[3] + 13 * maxnpair, Psize * Ssize, 0.0);
                                        std::fill_n(pTemps[3] + 14 * maxnpair, Psize * Ssize, 0.0);
                                        std::fill_n(pTemps[3] + 15 * maxnpair, Qsize * Rsize, 0.0);
                                        std::fill_n(pTemps[3] + 16 * maxnpair, Qsize * Rsize, 0.0);
                                        std::fill_n(pTemps[3] + 17 * maxnpair, Qsize * Rsize, 0.0);
                                    }
                                    delta = 0L;
                                    prefactor *= -0.25 * Kscale;
                                    for (int p = 0; p < Psize; p++) {
                                        for (int q = 0; q < Qsize; q++) {
                                            for (int r = 0; r < Rsize; r++) {
                                                for (int s = 0; s < Ssize; s++) {
                                                    Ax = prefactor * pAx[delta];
                                                    Ay = prefactor * pAy[delta];
                                                    Az = prefactor * pAz[delta];
                                                    Bx = prefactor * pBx[delta];
                                                    By = prefactor * pBy[delta];
                                                    Bz = prefactor * pBz[delta];
                                                    Cx = prefactor * pCx[delta];
                                                    Cy = prefactor * pCy[delta];
                                                    Cz = prefactor * pCz[delta];
                                                    Dx = prefactor * pDx[delta];
                                                    Dy = prefactor * pDy[delta];
                                                    Dz = prefactor * pDz[delta];

                                                    Dpr = Dap[p + Poff][r + Roff];
                                                    Dqs = Dap[q + Qoff][s + Soff];
                                                    Dps = Dap[p + Poff][s + Soff];
                                                    Dqr = Dap[q + Qoff][r + Roff];
                                                    if (pert_incore[Pcenter]) {
                                                        pTemps[0][6 * maxnpair + p * Rsize + r] += Ax * Dqs;
                                                        pTemps[0][7 * maxnpair + p * Rsize + r] += Ay * Dqs;
                                                        pTemps[0][8 * maxnpair + p * Rsize + r] += Az * Dqs;
                                                        pTemps[0][9 * maxnpair + q * Ssize + s] += Ax * Dpr;
                                                        pTemps[0][10 * maxnpair + q * Ssize + s] += Ay * Dpr;
                                                        pTemps[0][11 * maxnpair + q * Ssize + s] += Az * Dpr;
                                                        pTemps[0][12 * maxnpair + p * Ssize + s] += Ax * Dqr;
                                                        pTemps[0][13 * maxnpair + p * Ssize + s] += Ay * Dqr;
                                                        pTemps[0][14 * maxnpair + p * Ssize + s] += Az * Dqr;
                                                        pTemps[0][15 * maxnpair + q * Rsize + r] += Ax * Dps;
                                                        pTemps[0][16 * maxnpair + q * Rsize + r] += Ay * Dps;
                                                        pTemps[0][17 * maxnpair + q * Rsize + r] += Az * Dps;
                                                    }
                                                    if (pert_incore[Qcenter]) {
                                                        pTemps[1][6 * maxnpair + p * Rsize + r] += Bx * Dqs;
                                                        pTemps[1][7 * maxnpair + p * Rsize + r] += By * Dqs;
                                                        pTemps[1][8 * maxnpair + p * Rsize + r] += Bz * Dqs;
                                                        pTemps[1][9 * maxnpair + q * Ssize + s] += Bx * Dpr;
                                                        pTemps[1][10 * maxnpair + q * Ssize + s] += By * Dpr;
                                                        pTemps[1][11 * maxnpair + q * Ssize + s] += Bz * Dpr;
                                                        pTemps[1][12 * maxnpair + p * Ssize + s] += Bx * Dqr;
                                                        pTemps[1][13 * maxnpair + p * Ssize + s] += By * Dqr;
                                                        pTemps[1][14 * maxnpair + p * Ssize + s] += Bz * Dqr;
                                                        pTemps[1][15 * maxnpair + q * Rsize + r] += Bx * Dps;
                                                        pTemps[1][16 * maxnpair + q * Rsize + r] += By * Dps;
                                                        pTemps[1][17 * maxnpair + q * Rsize + r] += Bz * Dps;
                                                    }
                                                    if (pert_incore[Rcenter]) {
                                                        pTemps[2][6 * maxnpair + p * Rsize + r] += Cx * Dqs;
                                                        pTemps[2][7 * maxnpair + p * Rsize + r] += Cy * Dqs;
                                                        pTemps[2][8 * maxnpair + p * Rsize + r] += Cz * Dqs;
                                                        pTemps[2][9 * maxnpair + q * Ssize + s] += Cx * Dpr;
                                                        pTemps[2][10 * maxnpair + q * Ssize + s] += Cy * Dpr;
                                                        pTemps[2][11 * maxnpair + q * Ssize + s] += Cz * Dpr;
                                                        pTemps[2][12 * maxnpair + p * Ssize + s] += Cx * Dqr;
                                                        pTemps[2][13 * maxnpair + p * Ssize + s] += Cy * Dqr;
                                                        pTemps[2][14 * maxnpair + p * Ssize + s] += Cz * Dqr;
                                                        pTemps[2][15 * maxnpair + q * Rsize + r] += Cx * Dps;
                                                        pTemps[2][16 * maxnpair + q * Rsize + r] += Cy * Dps;
                                                        pTemps[2][17 * maxnpair + q * Rsize + r] += Cz * Dps;
                                                    }
                                                    if (pert_incore[Scenter]) {
                                                        pTemps[3][6 * maxnpair + p * Rsize + r] += Dx * Dqs;
                                                        pTemps[3][7 * maxnpair + p * Rsize + r] += Dy * Dqs;
                                                        pTemps[3][8 * maxnpair + p * Rsize + r] += Dz * Dqs;
                                                        pTemps[3][9 * maxnpair + q * Ssize + s] += Dx * Dpr;
                                                        pTemps[3][10 * maxnpair + q * Ssize + s] += Dy * Dpr;
                                                        pTemps[3][11 * maxnpair + q * Ssize + s] += Dz * Dpr;
                                                        pTemps[3][12 * maxnpair + p * Ssize + s] += Dx * Dqr;
                                                        pTemps[3][13 * maxnpair + p * Ssize + s] += Dy * Dqr;
                                                        pTemps[3][14 * maxnpair + p * Ssize + s] += Dz * Dqr;
                                                        pTemps[3][15 * maxnpair + q * Rsize + r] += Dx * Dps;
                                                        pTemps[3][16 * maxnpair + q * Rsize + r] += Dy * Dps;
                                                        pTemps[3][17 * maxnpair + q * Rsize + r] += Dz * Dps;
                                                    }
                                                    delta++;
                                                }
                                            }
                                        }
                                    }
                                }
                                // copy over temp data from per-thread buffers to the derivative fock matrices
                                if (pert_incore[Pcenter]) {
                                    // Jpq
                                    for (int xyz = 0; xyz < 3; ++xyz) {
                                        mutexes[3 * Pcenter + xyz].lock();
                                        for (int p = 0; p < Psize; ++p)
                                            for (int q = 0; q < Qsize; ++q)
                                                pdG[Pcenter * 3 + xyz][p + Poff][q + Qoff] +=
                                                    pTemps[0][(0 + xyz) * maxnpair + p * Qsize + q];
                                        mutexes[3 * Pcenter + xyz].unlock();
                                    }
                                    // Jrs
                                    for (int xyz = 0; xyz < 3; ++xyz) {
                                        mutexes[3 * Pcenter + xyz].lock();
                                        for (int r = 0; r < Rsize; ++r)
                                            for (int s = 0; s < Ssize; ++s)
                                                pdG[Pcenter * 3 + xyz][r + Roff][s + Soff] +=
                                                    pTemps[0][(3 + xyz) * maxnpair + r * Ssize + s];
                                        mutexes[3 * Pcenter + xyz].unlock();
                                    }
                                    if (Kscale != 0.0) {
                                        // Kpr
                                        for (int xyz = 0; xyz < 3; ++xyz) {
                                            mutexes[3 * Pcenter + xyz].lock();
                                            for (int p = 0; p < Psize; ++p)
                                                for (int r = 0; r < Rsize; ++r)
                                                    pdG[Pcenter * 3 + xyz][p + Poff][r + Roff] +=
                                                        pTemps[0][(6 + xyz) * maxnpair + p * Rsize + r];
                                            mutexes[3 * Pcenter + xyz].unlock();
                                        }
                                        // Kqs
                                        for (int xyz = 0; xyz < 3; ++xyz) {
                                            mutexes[3 * Pcenter + xyz].lock();
                                            for (int q = 0; q < Qsize; ++q)
                                                for (int s = 0; s < Ssize; ++s)
                                                    pdG[Pcenter * 3 + xyz][q + Qoff][s + Soff] +=
                                                        pTemps[0][(9 + xyz) * maxnpair + q * Ssize + s];
                                            mutexes[3 * Pcenter + xyz].unlock();
                                        }
                                        // Kps
                                        for (int xyz = 0; xyz < 3; ++xyz) {
                                            mutexes[3 * Pcenter + xyz].lock();
                                            for (int p = 0; p < Psize; ++p)
                                                for (int s = 0; s < Ssize; ++s)
                                                    pdG[Pcenter * 3 + xyz][p + Poff][s + Soff] +=
                                                        pTemps[0][(12 + xyz) * maxnpair + p * Ssize + s];
                                            mutexes[3 * Pcenter + xyz].unlock();
                                        }
                                        // Kqr
                                        for (int xyz = 0; xyz < 3; ++xyz) {
                                            mutexes[3 * Pcenter + xyz].lock();
                                            for (int q = 0; q < Qsize; ++q)
                                                for (int r = 0; r < Rsize; ++r)
                                                    pdG[Pcenter * 3 + xyz][q + Qoff][r + Roff] +=
                                                        pTemps[0][(15 + xyz) * maxnpair + q * Rsize + r];
                                            mutexes[3 * Pcenter + xyz].unlock();
                                        }
                                    }
                                }

                                if (pert_incore[Qcenter]) {
                                    // Jpq
                                    for (int xyz = 0; xyz < 3; ++xyz) {
                                        mutexes[3 * Qcenter + xyz].lock();
                                        for (int p = 0; p < Psize; ++p)
                                            for (int q = 0; q < Qsize; ++q)
                                                pdG[Qcenter * 3 + xyz][p + Poff][q + Qoff] +=
                                                    pTemps[1][(0 + xyz) * maxnpair + p * Qsize + q];
                                        mutexes[3 * Qcenter + xyz].unlock();
                                    }
                                    // Jrs
                                    for (int xyz = 0; xyz < 3; ++xyz) {
                                        mutexes[3 * Qcenter + xyz].lock();
                                        for (int r = 0; r < Rsize; ++r)
                                            for (int s = 0; s < Ssize; ++s)
                                                pdG[Qcenter * 3 + xyz][r + Roff][s + Soff] +=
                                                    pTemps[1][(3 + xyz) * maxnpair + r * Ssize + s];
                                        mutexes[3 * Qcenter + xyz].unlock();
                                    }
                                    if (Kscale != 0.0) {
                                        // Kpr
                                        for (int xyz = 0; xyz < 3; ++xyz) {
                                            mutexes[3 * Qcenter + xyz].lock();
                                            for (int p = 0; p < Psize; ++p)
                                                for (int r = 0; r < Rsize; ++r)
                                                    pdG[Qcenter * 3 + xyz][p + Poff][r + Roff] +=
                                                        pTemps[1][(6 + xyz) * maxnpair + p * Rsize + r];
                                            mutexes[3 * Qcenter + xyz].unlock();
                                        }
                                        // Kqs
                                        for (int xyz = 0; xyz < 3; ++xyz) {
                                            mutexes[3 * Qcenter + xyz].lock();
                                            for (int q = 0; q < Qsize; ++q)
                                                for (int s = 0; s < Ssize; ++s)
                                                    pdG[Qcenter * 3 + xyz][q + Qoff][s + Soff] +=
                                                        pTemps[1][(9 + xyz) * maxnpair + q * Ssize + s];
                                            mutexes[3 * Qcenter + xyz].unlock();
                                        }
                                        // Kps
                                        for (int xyz = 0; xyz < 3; ++xyz) {
                                            mutexes[3 * Qcenter + xyz].lock();
                                            for (int p = 0; p < Psize; ++p)
                                                for (int s = 0; s < Ssize; ++s)
                                                    pdG[Qcenter * 3 + xyz][p + Poff][s + Soff] +=
                                                        pTemps[1][(12 + xyz) * maxnpair + p * Ssize + s];
                                            mutexes[3 * Qcenter + xyz].unlock();
                                        }
                                        // Kqr
                                        for (int xyz = 0; xyz < 3; ++xyz) {
                                            mutexes[3 * Qcenter + xyz].lock();
                                            for (int q = 0; q < Qsize; ++q)
                                                for (int r = 0; r < Rsize; ++r)
                                                    pdG[Qcenter * 3 + xyz][q + Qoff][r + Roff] +=
                                                        pTemps[1][(15 + xyz) * maxnpair + q * Rsize + r];
                                            mutexes[3 * Qcenter + xyz].unlock();
                                        }
                                    }
                                }

                                if (pert_incore[Rcenter]) {
                                    // Jpq
                                    for (int xyz = 0; xyz < 3; ++xyz) {
                                        mutexes[3 * Rcenter + xyz].lock();
                                        for (int p = 0; p < Psize; ++p)
                                            for (int q = 0; q < Qsize; ++q)
                                                pdG[Rcenter * 3 + xyz][p + Poff][q + Qoff] +=
                                                    pTemps[2][(0 + xyz) * maxnpair + p * Qsize + q];
                                        mutexes[3 * Rcenter + xyz].unlock();
                                    }
                                    // Jrs
                                    for (int xyz = 0; xyz < 3; ++xyz) {
                                        mutexes[3 * Rcenter + xyz].lock();
                                        for (int r = 0; r < Rsize; ++r)
                                            for (int s = 0; s < Ssize; ++s)
                                                pdG[Rcenter * 3 + xyz][r + Roff][s + Soff] +=
                                                    pTemps[2][(3 + xyz) * maxnpair + r * Ssize + s];
                                        mutexes[3 * Rcenter + xyz].unlock();
                                    }
                                    if (Kscale != 0.0) {
                                        // Kpr
                                        for (int xyz = 0; xyz < 3; ++xyz) {
                                            mutexes[3 * Rcenter + xyz].lock();
                                            for (int p = 0; p < Psize; ++p)
                                                for (int r = 0; r < Rsize; ++r)
                                                    pdG[Rcenter * 3 + xyz][p + Poff][r + Roff] +=
                                                        pTemps[2][(6 + xyz) * maxnpair + p * Rsize + r];
                                            mutexes[3 * Rcenter + xyz].unlock();
                                        }
                                        // Kqs
                                        for (int xyz = 0; xyz < 3; ++xyz) {
                                            mutexes[3 * Rcenter + xyz].lock();
                                            for (int q = 0; q < Qsize; ++q)
                                                for (int s = 0; s < Ssize; ++s)
                                                    pdG[Rcenter * 3 + xyz][q + Qoff][s + Soff] +=
                                                        pTemps[2][(9 + xyz) * maxnpair + q * Ssize + s];
                                            mutexes[3 * Rcenter + xyz].unlock();
                                        }
                                        // Kps
                                        for (int xyz = 0; xyz < 3; ++xyz) {
                                            mutexes[3 * Rcenter + xyz].lock();
                                            for (int p = 0; p < Psize; ++p)
                                                for (int s = 0; s < Ssize; ++s)
                                                    pdG[Rcenter * 3 + xyz][p + Poff][s + Soff] +=
                                                        pTemps[2][(12 + xyz) * maxnpair + p * Ssize + s];
                                            mutexes[3 * Rcenter + xyz].unlock();
                                        }
                                        // Kqr
                                        for (int xyz = 0; xyz < 3; ++xyz) {
                                            mutexes[3 * Rcenter + xyz].lock();
                                            for (int q = 0; q < Qsize; ++q)
                                                for (int r = 0; r < Rsize; ++r)
                                                    pdG[Rcenter * 3 + xyz][q + Qoff][r + Roff] +=
                                                        pTemps[2][(15 + xyz) * maxnpair + q * Rsize + r];
                                            mutexes[3 * Rcenter + xyz].unlock();
                                        }
                                    }
                                }
                                if (pert_incore[Scenter]) {
                                    // Jpq
                                    for (int xyz = 0; xyz < 3; ++xyz) {
                                        mutexes[3 * Scenter + xyz].lock();
                                        for (int p = 0; p < Psize; ++p)
                                            for (int q = 0; q < Qsize; ++q)
                                                pdG[Scenter * 3 + xyz][p + Poff][q + Qoff] +=
                                                    pTemps[3][(0 + xyz) * maxnpair + p * Qsize + q];
                                        mutexes[3 * Scenter + xyz].unlock();
                                    }
                                    // Jrs
                                    for (int xyz = 0; xyz < 3; ++xyz) {
                                        mutexes[3 * Scenter + xyz].lock();
                                        for (int r = 0; r < Rsize; ++r)
                                            for (int s = 0; s < Ssize; ++s)
                                                pdG[Scenter * 3 + xyz][r + Roff][s + Soff] +=
                                                    pTemps[3][(3 + xyz) * maxnpair + r * Ssize + s];
                                        mutexes[3 * Scenter + xyz].unlock();
                                    }
                                    if (Kscale != 0.0) {
                                        // Kpr
                                        for (int xyz = 0; xyz < 3; ++xyz) {
                                            mutexes[3 * Scenter + xyz].lock();
                                            for (int p = 0; p < Psize; ++p)
                                                for (int r = 0; r < Rsize; ++r)
                                                    pdG[Scenter * 3 + xyz][p + Poff][r + Roff] +=
                                                        pTemps[3][(6 + xyz) * maxnpair + p * Rsize + r];
                                            mutexes[3 * Scenter + xyz].unlock();
                                        }
                                        // Kqs
                                        for (int xyz = 0; xyz < 3; ++xyz) {
                                            mutexes[3 * Scenter + xyz].lock();
                                            for (int q = 0; q < Qsize; ++q)
                                                for (int s = 0; s < Ssize; ++s)
                                                    pdG[Scenter * 3 + xyz][q + Qoff][s + Soff] +=
                                                        pTemps[3][(9 + xyz) * maxnpair + q * Ssize + s];
                                            mutexes[3 * Scenter + xyz].unlock();
                                        }
                                        // Kps
                                        for (int xyz = 0; xyz < 3; ++xyz) {
                                            mutexes[3 * Scenter + xyz].lock();
                                            for (int p = 0; p < Psize; ++p)
                                                for (int s = 0; s < Ssize; ++s)
                                                    pdG[Scenter * 3 + xyz][p + Poff][s + Soff] +=
                                                        pTemps[3][(12 + xyz) * maxnpair + p * Ssize + s];
                                            mutexes[3 * Scenter + xyz].unlock();
                                        }
                                        // Kqr
                                        for (int xyz = 0; xyz < 3; ++xyz) {
                                            mutexes[3 * Scenter + xyz].lock();
                                            for (int q = 0; q < Qsize; ++q)
                                                for (int r = 0; r < Rsize; ++r)
                                                    pdG[Scenter * 3 + xyz][q + Qoff][r + Roff] +=
                                                        pTemps[3][(15 + xyz) * maxnpair + q * Rsize + r];
                                            mutexes[3 * Scenter + xyz].unlock();
                                        }
                                    }
                                }

                                pAx += block_size;
                                pAy += block_size;
                                pAz += block_size;
                                pBx += block_size;
                                pBy += block_size;
                                pBz += block_size;
                                pCx += block_size;
                                pCy += block_size;
                                pCz += block_size;
                                pDx += block_size;
                                pDy += block_size;
                                pDz += block_size;
                            }  // pairRS
                        }      // pairPQ
                    }          // blockRS
                }              // blockPQ

                for (int a = 0; a < nA; ++a) {
                    // Symmetrize the derivative Fock contributions
                    SharedMatrix G = dGmats[a];
                    G->add(G->transpose());
                    Gpi->transform(C, G, Cocc);
                    Gpi->scale(0.5);
                    psio_->write(PSIF_HESS, "Gpi^A", (char*)pGpi[0], static_cast<size_t>(nmo) * nocc * sizeof(double),
                                 next_Gpi, &next_Gpi);
                }
            } // End loop over A batches
#ifdef USING_BrianQC
            }
#endif
        } // End if density fitted
    }

    size_t mem = 0.9 * memory_ / 8L;
    size_t per_A = 3L * nso * nso + 1L * nocc * nso;
    size_t max_A = (mem / 2L) / per_A;
    max_A = (max_A > 3 * natom ? 3 * natom : max_A);

    std::shared_ptr<JK> jk;
    jk = JK::build_JK(basisset_, get_basisset("DF_BASIS_SCF"), options_, false, mem);

    jk->set_memory(mem);
    jk->initialize();

    // => J2pi/K2pi <= //
    {
        // Figure out DFT functional info
        double Kscale = functional_->x_alpha();
        if (functional_->is_x_lrc()) throw PSIEXCEPTION("Hessians for LRC functionals are not implemented yet.");

        std::vector<std::shared_ptr<Matrix>>& L = jk->C_left();
        std::vector<std::shared_ptr<Matrix>>& R = jk->C_right();
        const std::vector<std::shared_ptr<Matrix>>& J = jk->J();
        const std::vector<std::shared_ptr<Matrix>>& K = jk->K();
        L.clear();
        R.clear();

        auto Sij = std::make_shared<Matrix>("Sij", nocc, nocc);
        double** Sijp = Sij->pointer();
        auto T = std::make_shared<Matrix>("T", nso, nocc);
        double** Tp = T->pointer();
        auto U = std::make_shared<Matrix>("Tempai", nmo, nocc);
        double** Up = U->pointer();

        // Write some placeholder data to PSIO, to get the sizing right
        psio_address next_Gpi = PSIO_ZERO;
        for (int A = 0; A < 3 * natom; A++)
            psio_->write(PSIF_HESS, "G2pi^A", (char*)Up[0], static_cast<size_t>(nmo) * nocc * sizeof(double), next_Gpi,
                         &next_Gpi);

        std::vector<SharedMatrix> Dx, Vx;
        for (int A = 0; A < max_A; A++) {
            // Just pass C1 quantities in; this object doesn't respect symmetry anyway
            L.push_back(Cocc);
            R.push_back(std::make_shared<Matrix>("R", nso, nocc));
            Dx.push_back(std::make_shared<Matrix>("Dx", nso, nso));
            Vx.push_back(std::make_shared<Matrix>("Vx", nso, nso));
        }

        jk->print_header();

        for (int A = 0; A < 3 * natom; A += max_A) {
            int nA = max_A;
            if (A + max_A >= 3 * natom) {
                nA = 3 * natom - A;
                L.resize(nA);
                R.resize(nA);
            }
            for (int a = 0; a < nA; a++) {
                psio_address next_Sij = psio_get_address(PSIO_ZERO, (A + a) * (size_t)nocc * nocc * sizeof(double));
                psio_->read(PSIF_HESS, "Sij^A", (char*)Sijp[0], static_cast<size_t>(nocc) * nocc * sizeof(double),
                            next_Sij, &next_Sij);
                C_DGEMM('N', 'N', nso, nocc, nocc, 1.0, Cop[0], nocc, Sijp[0], nocc, 0.0, R[a]->pointer()[0], nocc);
                Dx[a] = linalg::doublet(L[a], R[a], false, true);
                // Symmetrize the pseudodensity
                Dx[a]->add(Dx[a]->transpose());
                Dx[a]->scale(0.5);
            }

            jk->compute();
            if (functional_->needs_xc()) {
                potential_->compute_Vx(Dx, Vx);
            }

            for (int a = 0; a < nA; a++) {
                // Add the 2J contribution to G
                C_DGEMM('N', 'N', nso, nocc, nso, 1.0, J[a]->pointer()[0], nso, Cop[0], nocc, 0.0, Tp[0], nocc);
                C_DGEMM('T', 'N', nmo, nocc, nso, -2.0, Cp[0], nmo, Tp[0], nocc, 0.0, Up[0], nocc);

                if (functional_->needs_xc()) {
                    // Symmetrize the result, just to be safe
                    C_DGEMM('N', 'N', nso, nocc, nso, 0.5, Vx[a]->pointer()[0], nso, Cop[0], nocc, 0.0, Tp[0], nocc);
                    C_DGEMM('T', 'N', nso, nocc, nso, 0.5, Vx[a]->pointer()[0], nso, Cop[0], nocc, 1.0, Tp[0], nocc);
                    C_DGEMM('T', 'N', nmo, nocc, nso, -2.0, Cp[0], nmo, Tp[0], nocc, 1.0, Up[0], nocc);
                }

                // Subtract the K term from G
                if (Kscale) {
                    C_DGEMM('N', 'N', nso, nocc, nso, 1.0, K[a]->pointer()[0], nso, Cop[0], nocc, 0.0, Tp[0], nocc);
                    C_DGEMM('T', 'N', nmo, nocc, nso, Kscale, Cp[0], nmo, Tp[0], nocc, 1.0, Up[0], nocc);
                }

                psio_address next_Gpi = psio_get_address(PSIO_ZERO, (A + a) * (size_t)nmo * nocc * sizeof(double));
                psio_->write(PSIF_HESS, "G2pi^A", (char*)Up[0], static_cast<size_t>(nmo) * nocc * sizeof(double),
                             next_Gpi, &next_Gpi);
            }
        }
    }

    // => XC Gradient <= //
    {
        auto T = std::make_shared<Matrix>("T", nso, nocc);
        double** Tp = T->pointer();
        auto U = std::make_shared<Matrix>("Tempai", nmo, nocc);
        double** Up = U->pointer();
        if (functional_->needs_xc()) {
            // Write some placeholder data to PSIO, to get the sizing right
            psio_address next_VXCpi = PSIO_ZERO;
            for (int A = 0; A < 3 * natom; A++)
                psio_->write(PSIF_HESS, "VXCpi^A", (char*)Up[0], static_cast<size_t>(nmo) * nocc * sizeof(double),
                             next_VXCpi, &next_VXCpi);
            // For now we just compute all 3N matrices in one go.  If this becomes to burdensome
            // in terms of memory we can reevaluate and implement a batching mechanism instead.
            auto Vxc_matrices = potential_->compute_fock_derivatives();
            for (int a = 0; a < 3 * natom; ++a) {
                // Transform from SO basis to pi
                C_DGEMM('N', 'N', nso, nocc, nso, 1.0, Vxc_matrices[a]->pointer()[0], nso, Cop[0], nocc, 0.0, Tp[0],
                        nocc);
                C_DGEMM('T', 'N', nmo, nocc, nso, 1.0, Cp[0], nmo, Tp[0], nocc, 0.0, Up[0], nocc);
                next_VXCpi = psio_get_address(PSIO_ZERO, a * (size_t)nmo * nocc * sizeof(double));
                psio_->write(PSIF_HESS, "VXCpi^A", (char*)Up[0], static_cast<size_t>(nmo) * nocc * sizeof(double),
                             next_VXCpi, &next_VXCpi);
            }
        }
    }

    // => Fpi <= //
    {
        auto Tpi = std::make_shared<Matrix>("Tpi", nmo, nocc);
        auto Fpi = std::make_shared<Matrix>("Fpi", nmo, nocc);
        double** Tpip = Tpi->pointer();
        double** Fpip = Fpi->pointer();

        psio_address next_Epi = PSIO_ZERO;
        psio_address next_Tpi = PSIO_ZERO;
        psio_address next_Vpi = PSIO_ZERO;
        psio_address next_Jpi = PSIO_ZERO;
        psio_address next_Fpi = PSIO_ZERO;
        psio_address next_VXCpi = PSIO_ZERO;

        for (int A = 0; A < 3 * natom; A++) {
            psio_->read(PSIF_HESS, "Tpi^A", (char*)Fpip[0], static_cast<size_t>(nmo) * nocc * sizeof(double), next_Tpi,
                        &next_Tpi);
            psio_->read(PSIF_HESS, "Vpi^A", (char*)Tpip[0], static_cast<size_t>(nmo) * nocc * sizeof(double), next_Vpi,
                        &next_Vpi);
            Fpi->add(Tpi);
            if (basisset_->has_ECP()) {
                psio_->read(PSIF_HESS, "Epi^A", (char*)Tpip[0], static_cast<size_t>(nmo) * nocc * sizeof(double), next_Epi,
                            &next_Epi);
                Fpi->add(Tpi);
            }
            psio_->read(PSIF_HESS, "Gpi^A", (char*)Tpip[0], static_cast<size_t>(nmo) * nocc * sizeof(double), next_Jpi,
                        &next_Jpi);
            Fpi->add(Tpi);
            if (functional_->needs_xc()) {
                psio_->read(PSIF_HESS, "VXCpi^A", (char*)Tpip[0], static_cast<size_t>(nmo) * nocc * sizeof(double),
                            next_VXCpi, &next_VXCpi);
                Fpi->add(Tpi);
            }
            psio_->write(PSIF_HESS, "Fpi^A", (char*)Fpip[0], static_cast<size_t>(nmo) * nocc * sizeof(double), next_Fpi,
                         &next_Fpi);
        }
    }

    // => Bpi <= //
    {
        auto Tai = std::make_shared<Matrix>("T", nvir, nocc);
        auto Bai = std::make_shared<Matrix>("B", nvir, nocc);
        double** Taip = Tai->pointer();
        double** Baip = Bai->pointer();

        psio_address next_Fpi = PSIO_ZERO;
        psio_address next_Spi = PSIO_ZERO;
        psio_address next_G2pi = PSIO_ZERO;
        psio_address next_Bai = PSIO_ZERO;

        for (int A = 0; A < 3 * natom; A++) {
            next_Fpi = psio_get_address(PSIO_ZERO, sizeof(double) * (A * (size_t)nmo * nocc + nocc * nocc));
            psio_->read(PSIF_HESS, "Fpi^A", (char*)Baip[0], static_cast<size_t>(nvir) * nocc * sizeof(double), next_Fpi,
                        &next_Fpi);
            next_Spi = psio_get_address(PSIO_ZERO, sizeof(double) * (A * (size_t)nmo * nocc + nocc * nocc));
            psio_->read(PSIF_HESS, "Spi^A", (char*)Taip[0], static_cast<size_t>(nvir) * nocc * sizeof(double), next_Spi,
                        &next_Spi);
            for (int i = 0; i < nocc; i++) C_DAXPY(nvir, -eop[i], &Taip[0][i], nocc, &Baip[0][i], nocc);
            next_G2pi = psio_get_address(PSIO_ZERO, sizeof(double) * (A * (size_t)nmo * nocc + nocc * nocc));
            psio_->read(PSIF_HESS, "G2pi^A", (char*)Taip[0], static_cast<size_t>(nvir) * nocc * sizeof(double),
                        next_G2pi, &next_G2pi);
            Bai->add(Tai);
            Bai->scale(-1.0);
            psio_->write(PSIF_HESS, "Bai^A", (char*)Baip[0], static_cast<size_t>(nvir) * nocc * sizeof(double),
                         next_Bai, &next_Bai);
        }
    }

    // => CPHF (Uai) <= //
    {
        rhf_wfn_->set_jk(jk);

        psio_address next_Bai = PSIO_ZERO;
        psio_address next_Uai = PSIO_ZERO;

        auto T = std::make_shared<Matrix>("T", nvir, nocc);
        double** Tp = T->pointer();

        for (int A = 0; A < 3 * natom; A += max_A) {
            int nA = max_A;
            if (A + max_A >= 3 * natom) {
                nA = 3 * natom - A;
            }

            std::vector<SharedMatrix> b_vecs;
            // Fill b
            for (int a = 0; a < nA; a++) {
                std::stringstream ss;
                ss << "Perturbation " << a + A;
                auto B = std::make_shared<Matrix>(ss.str(), nocc, nvir);
                psio_->read(PSIF_HESS, "Bai^A", (char*)Tp[0], static_cast<size_t>(nvir) * nocc * sizeof(double),
                            next_Bai, &next_Bai);
                double** Bp = B->pointer();
                for (int i = 0; i < nocc; i++) {
                    C_DCOPY(nvir, &Tp[0][i], nocc, Bp[i], 1);
                }
                b_vecs.push_back(B);
            }

#ifdef USING_BrianQC
    brianCPHFLeftSideFlag = true;
#endif
            auto u_matrices = rhf_wfn_->cphf_solve(b_vecs, options_.get_double("SOLVER_CONVERGENCE"),
                                                   options_.get_int("SOLVER_MAXITER"), print_);

#ifdef USING_BrianQC
    brianCPHFLeftSideFlag = false;
#endif
            // Result in x
            for (int a = 0; a < nA; a++) {
                std::stringstream ss;
                ss << "Perturbation " << a + A;
                u_matrices[a]->scale(-1);
                double** Xp = u_matrices[a]->pointer();
                for (int i = 0; i < nocc; i++) {
                    C_DCOPY(nvir, Xp[i], 1, &Tp[0][i], nocc);
                }
                psio_->write(PSIF_HESS, "Uai^A", (char*)Tp[0], static_cast<size_t>(nvir) * nocc * sizeof(double),
                             next_Uai, &next_Uai);
            }
        }
    }

    // => Dipole derivatives (for IR intensities) <= //
    MintsHelper mints(basisset_);
    auto ao_dipole = mints.ao_dipole();
    auto Ca = Ca_subset("AO");
    auto Caocc = Ca_subset("AO", "OCC");
    Matrix mu_x("mu X", nmo_, nocc);
    Matrix mu_y("mu Y", nmo_, nocc);
    Matrix mu_z("mu Z", nmo_, nocc);
    mu_x.transform(Ca, ao_dipole[0], Caocc);
    mu_y.transform(Ca, ao_dipole[1], Caocc);
    mu_z.transform(Ca, ao_dipole[2], Caocc);
    // Start by computing the skeleton derivatives
    auto dipole_gradient = mints.dipole_grad(Da_subset("AO"));
    // Account for alpha and beta orbitals
    dipole_gradient->scale(2);
    dipole_gradient->add(DipoleInt::nuclear_gradient_contribution(molecule_));
    double** pdip_grad = dipole_gradient->pointer();

    // => Upi <= //
    {
        auto Upi = std::make_shared<Matrix>("U", nmo, nocc);
        double** Upqp = Upi->pointer();

        psio_address next_Spi = PSIO_ZERO;
        psio_address next_Uai = PSIO_ZERO;
        psio_address next_Upi = PSIO_ZERO;

        for (int A = 0; A < 3 * natom; A++) {
            psio_->read(PSIF_HESS, "Sij^A", (char*)Upqp[0], static_cast<size_t>(nocc) * nocc * sizeof(double), next_Spi,
                        &next_Spi);
            C_DSCAL(nocc * (size_t)nocc, -0.5, Upqp[0], 1);
            psio_->read(PSIF_HESS, "Uai^A", (char*)Upqp[nocc], static_cast<size_t>(nvir) * nocc * sizeof(double),
                        next_Uai, &next_Uai);
            psio_->write(PSIF_HESS, "Upi^A", (char*)Upqp[0], static_cast<size_t>(nmo) * nocc * sizeof(double), next_Upi,
                         &next_Upi);
            pdip_grad[A][0] += 4 * mu_x.vector_dot(Upi);
            pdip_grad[A][1] += 4 * mu_y.vector_dot(Upi);
            pdip_grad[A][2] += 4 * mu_z.vector_dot(Upi);
        }
    }
    rhf_wfn_->set_array_variable("SCF DIPOLE GRADIENT", dipole_gradient);
    rhf_wfn_->set_array_variable("CURRENT DIPOLE GRADIENT", dipole_gradient);

    // => Qpi <= //
    {
        double Kscale = functional_->x_alpha();
        std::vector<std::shared_ptr<Matrix>>& L = jk->C_left();
        std::vector<std::shared_ptr<Matrix>>& R = jk->C_right();
        const std::vector<std::shared_ptr<Matrix>>& J = jk->J();
        const std::vector<std::shared_ptr<Matrix>>& K = jk->K();
        std::vector<SharedMatrix> Dx, Vx;
        L.clear();
        R.clear();
        for (int a = 0; a < max_A; a++) {
            L.push_back(Cocc);
            R.push_back(std::make_shared<Matrix>("R", nso, nocc));
            Dx.push_back(std::make_shared<Matrix>("Dx", nso, nso));
            Vx.push_back(std::make_shared<Matrix>("Vx", nso, nso));
        }

        auto Upi = std::make_shared<Matrix>("Upi", nmo, nocc);
        double** Upip = Upi->pointer();
        auto T = std::make_shared<Matrix>("T", nso, nocc);
        double** Tp = T->pointer();
        auto U = std::make_shared<Matrix>("T", nmo, nocc);
        double** Up = U->pointer();

        for (int A = 0; A < 3 * natom; A += max_A) {
            int nA = max_A;
            if (A + max_A >= 3 * natom) {
                nA = 3 * natom - A;
                L.resize(nA);
                R.resize(nA);
            }
            for (int a = 0; a < nA; a++) {
                psio_address next_Upi = psio_get_address(PSIO_ZERO, (A + a) * (size_t)nmo * nocc * sizeof(double));
                psio_->read(PSIF_HESS, "Upi^A", (char*)Upip[0], static_cast<size_t>(nmo) * nocc * sizeof(double),
                            next_Upi, &next_Upi);
                C_DGEMM('N', 'N', nso, nocc, nmo, 1.0, Cp[0], nmo, Upip[0], nocc, 0.0, R[a]->pointer()[0], nocc);
                Dx[a] = linalg::doublet(L[a], R[a], false, true);
                // Symmetrize the pseudodensity
                Dx[a]->add(Dx[a]->transpose());
                Dx[a]->scale(0.5);
            }

            jk->compute();
            if (functional_->needs_xc()) {
                potential_->compute_Vx(Dx, Vx);
            }

            for (int a = 0; a < nA; a++) {
                C_DGEMM('N', 'N', nso, nocc, nso, 4.0, J[a]->pointer()[0], nso, Cop[0], nocc, 0.0, Tp[0], nocc);
                if (Kscale) {
                    C_DGEMM('N', 'N', nso, nocc, nso, -Kscale, K[a]->pointer()[0], nso, Cop[0], nocc, 1.0, Tp[0], nocc);
                    C_DGEMM('T', 'N', nso, nocc, nso, -Kscale, K[a]->pointer()[0], nso, Cop[0], nocc, 1.0, Tp[0], nocc);
                }
                if (functional_->needs_xc()) {
                    // Symmetrize the result, just to be safe
                    C_DGEMM('N', 'N', nso, nocc, nso, 2.0, Vx[a]->pointer()[0], nso, Cop[0], nocc, 1.0, Tp[0], nocc);
                    C_DGEMM('T', 'N', nso, nocc, nso, 2.0, Vx[a]->pointer()[0], nso, Cop[0], nocc, 1.0, Tp[0], nocc);
                }
                C_DGEMM('T', 'N', nmo, nocc, nso, 1.0, Cp[0], nmo, Tp[0], nocc, 0.0, Up[0], nocc);
                psio_address next_Qpi = psio_get_address(PSIO_ZERO, (A + a) * (size_t)nmo * nocc * sizeof(double));
                psio_->write(PSIF_HESS, "Qpi^A", (char*)Up[0], static_cast<size_t>(nmo) * nocc * sizeof(double),
                             next_Qpi, &next_Qpi);
            }
        }
    }
    jk.reset();

    // => Zipper <= //
    {
        size_t memory = 0.9 * memory_ / 8L;
        size_t npi = nmo * (size_t)nocc;
        size_t max_a = memory / (3L * npi);
        max_a = (max_a > 3 * natom ? 3 * natom : max_a);

        auto L = std::make_shared<Matrix>("L", max_a * nmo, nocc);
        auto R = std::make_shared<Matrix>("R", max_a * nmo, nocc);
        auto T = std::make_shared<Matrix>("T", max_a * nmo, nocc);
        double** Lp = L->pointer();
        double** Rp = R->pointer();
        double** Tp = T->pointer();

        double** Hp = response->pointer();

        // U^A F^B
        for (int A = 0; A < 3 * natom; A += max_a) {
            int nA = (A + max_a >= 3 * natom ? 3 * natom - A : max_a);
            psio_address nextA = psio_get_address(PSIO_ZERO, A * npi * sizeof(double));
            psio_->read(PSIF_HESS, "Upi^A", (char*)Lp[0], sizeof(double) * nA * npi, nextA, &nextA);
            for (int B = 0; B < 3 * natom; B += max_a) {
                int nB = (B + max_a >= 3 * natom ? 3 * natom - B : max_a);
                psio_address nextB = psio_get_address(PSIO_ZERO, B * npi * sizeof(double));
                psio_->read(PSIF_HESS, "Fpi^A", (char*)Rp[0], sizeof(double) * nB * npi, nextB, &nextB);
                for (int a = 0; a < nA; a++) {
                    for (int b = 0; b < nB; b++) {
                        Hp[A + a][B + b] += 4.0 * C_DDOT(npi, Lp[0] + a * npi, 1, Rp[0] + b * npi, 1);
                    }
                }
            }
        }

        // F^A U^B
        for (int A = 0; A < 3 * natom; A += max_a) {
            int nA = (A + max_a >= 3 * natom ? 3 * natom - A : max_a);
            psio_address nextA = psio_get_address(PSIO_ZERO, A * npi * sizeof(double));
            psio_->read(PSIF_HESS, "Fpi^A", (char*)Lp[0], sizeof(double) * nA * npi, nextA, &nextA);
            for (int B = 0; B < 3 * natom; B += max_a) {
                int nB = (B + max_a >= 3 * natom ? 3 * natom - B : max_a);
                psio_address nextB = psio_get_address(PSIO_ZERO, B * npi * sizeof(double));
                psio_->read(PSIF_HESS, "Upi^A", (char*)Rp[0], sizeof(double) * nB * npi, nextB, &nextB);
                for (int a = 0; a < nA; a++) {
                    for (int b = 0; b < nB; b++) {
                        Hp[A + a][B + b] += 4.0 * C_DDOT(npi, Lp[0] + a * npi, 1, Rp[0] + b * npi, 1);
                    }
                }
            }
        }

        // U^A U^B
        // N.B. We use the relationship U^a_ia = -U^a_ai - S^a_ai
        psio_address junk;
        for (int A = 0; A < 3 * natom; A += max_a) {
            int nA = (A + max_a >= 3 * natom ? 3 * natom - A : max_a);
            psio_address nextA = psio_get_address(PSIO_ZERO, A * npi * sizeof(double));
            psio_->read(PSIF_HESS, "Upi^A", (char*)Lp[0], sizeof(double) * nA * npi, nextA, &junk);
            psio_->read(PSIF_HESS, "Spi^A", (char*)Tp[0], sizeof(double) * nA * npi, nextA, &junk);
            L->add(T);
            for (int i = 0; i < nocc; i++) C_DSCAL(static_cast<size_t>(nA) * nmo, eop[i], &Lp[0][i], nocc);
            for (int B = 0; B < 3 * natom; B += max_a) {
                int nB = (B + max_a >= 3 * natom ? 3 * natom - B : max_a);
                psio_address nextB = psio_get_address(PSIO_ZERO, B * npi * sizeof(double));
                psio_->read(PSIF_HESS, "Upi^A", (char*)Rp[0], sizeof(double) * nB * npi, nextB, &junk);
                psio_->read(PSIF_HESS, "Spi^A", (char*)Tp[0], sizeof(double) * nB * npi, nextB, &junk);
                R->add(T);
                for (int a = 0; a < nA; a++) {
                    for (int b = 0; b < nB; b++) {
                        Hp[A + a][B + b] -= 2.0 * C_DDOT(npi, Lp[0] + a * npi, 1, Rp[0] + b * npi, 1);
                        Hp[B + b][A + a] -= 2.0 * C_DDOT(npi, Lp[0] + a * npi, 1, Rp[0] + b * npi, 1);
                    }
                }
            }
        }

        // S^A S^B
        for (int A = 0; A < 3 * natom; A += max_a) {
            int nA = (A + max_a >= 3 * natom ? 3 * natom - A : max_a);
            psio_address nextA = psio_get_address(PSIO_ZERO, A * npi * sizeof(double));
            psio_->read(PSIF_HESS, "Spi^A", (char*)Lp[0], sizeof(double) * nA * npi, nextA, &nextA);
            for (int i = 0; i < nocc; i++) C_DSCAL(static_cast<size_t>(nA) * nmo, eop[i], &Lp[0][i], nocc);
            for (int B = 0; B < 3 * natom; B += max_a) {
                int nB = (B + max_a >= 3 * natom ? 3 * natom - B : max_a);
                psio_address nextB = psio_get_address(PSIO_ZERO, B * npi * sizeof(double));
                psio_->read(PSIF_HESS, "Spi^A", (char*)Rp[0], sizeof(double) * nB * npi, nextB, &nextB);
                for (int a = 0; a < nA; a++) {
                    for (int b = 0; b < nB; b++) {
                        Hp[A + a][B + b] += 2.0 * C_DDOT(npi, Lp[0] + a * npi, 1, Rp[0] + b * npi, 1);
                        Hp[B + b][A + a] += 2.0 * C_DDOT(npi, Lp[0] + a * npi, 1, Rp[0] + b * npi, 1);
                    }
                }
            }
        }

        // U^A U^B \epsilon
        for (int A = 0; A < 3 * natom; A += max_a) {
            int nA = (A + max_a >= 3 * natom ? 3 * natom - A : max_a);
            psio_address nextA = psio_get_address(PSIO_ZERO, A * npi * sizeof(double));
            psio_->read(PSIF_HESS, "Upi^A", (char*)Lp[0], sizeof(double) * nA * npi, nextA, &nextA);
            double* Tp = Lp[0];
            for (int a = 0; a < nA; a++) {
                for (int p = 0; p < nmo; p++) {
                    C_DSCAL(nocc, ep[p], Tp, 1);
                    Tp += nocc;
                }
            }
            for (int B = 0; B < 3 * natom; B += max_a) {
                int nB = (B + max_a >= 3 * natom ? 3 * natom - B : max_a);
                psio_address nextB = psio_get_address(PSIO_ZERO, B * npi * sizeof(double));
                psio_->read(PSIF_HESS, "Upi^A", (char*)Rp[0], sizeof(double) * nB * npi, nextB, &nextB);
                for (int a = 0; a < nA; a++) {
                    for (int b = 0; b < nB; b++) {
                        Hp[A + a][B + b] += 4.0 * C_DDOT(npi, Lp[0] + a * npi, 1, Rp[0] + b * npi, 1);
                    }
                }
            }
        }

        // U^A Q^B
        for (int A = 0; A < 3 * natom; A += max_a) {
            int nA = (A + max_a >= 3 * natom ? 3 * natom - A : max_a);
            psio_address nextA = psio_get_address(PSIO_ZERO, A * npi * sizeof(double));
            psio_->read(PSIF_HESS, "Upi^A", (char*)Lp[0], sizeof(double) * nA * npi, nextA, &nextA);
            for (int B = 0; B < 3 * natom; B += max_a) {
                int nB = (B + max_a >= 3 * natom ? 3 * natom - B : max_a);
                psio_address nextB = psio_get_address(PSIO_ZERO, B * npi * sizeof(double));
                psio_->read(PSIF_HESS, "Qpi^A", (char*)Rp[0], sizeof(double) * nB * npi, nextB, &nextB);
                for (int a = 0; a < nA; a++) {
                    for (int b = 0; b < nB; b++) {
                        Hp[A + a][B + b] += 4.0 * C_DDOT(npi, Lp[0] + a * npi, 1, Rp[0] + b * npi, 1);
                    }
                }
            }
        }

        // Full symmetrization
        for (int A = 0; A < 3 * natom; A++) {
            for (int B = 0; B < 3 * natom; B++) {
                Hp[A][B] = Hp[B][A] = 0.5 * (Hp[A][B] + Hp[B][A]);
            }
        }
    }

    psio_->close(PSIF_HESS, 0);
    return response;
}

// Basically copy/paste the rhf function, split into spin compnents
std::shared_ptr<Matrix> USCFDeriv::hessian_response()
{
    // => Control Parameters <= //
    std::shared_ptr<Vector> eps_a     = epsilon_a_subset("AO","ALL");
    std::shared_ptr<Vector> eps_aocc = epsilon_a_subset("AO","OCC");
    std::shared_ptr<Vector> eps_avir = epsilon_a_subset("AO","VIR");
    std::shared_ptr<Matrix> Ca    = Ca_subset("AO","ALL");
    std::shared_ptr<Matrix> Ca_occ = Ca_subset("AO","OCC");
    std::shared_ptr<Matrix> Ca_vir = Ca_subset("AO","VIR");
    std::shared_ptr<Matrix> Da = Da_subset("AO");

    std::shared_ptr<Vector> eps_b     = epsilon_b_subset("AO","ALL");
    std::shared_ptr<Vector> eps_bocc = epsilon_b_subset("AO","OCC");
    std::shared_ptr<Vector> eps_bvir = epsilon_b_subset("AO","VIR");
    std::shared_ptr<Matrix> Cb    = Cb_subset("AO","ALL");
    std::shared_ptr<Matrix> Cb_occ = Cb_subset("AO","OCC");
    std::shared_ptr<Matrix> Cb_vir = Cb_subset("AO","VIR");
    std::shared_ptr<Matrix> Db = Db_subset("AO");

    std::shared_ptr<Matrix> Dt = Da_subset("AO")->clone();
    Dt->add(Db);

    // => Sizing <= //

    int natom = molecule_->natom();
    int nso   = basisset_->nbf();
    int naocc  = eps_aocc->dimpi()[0];
    int navir  = eps_avir->dimpi()[0];
    int nbocc  = eps_bocc->dimpi()[0];
    int nbvir  = eps_bvir->dimpi()[0];
    int nmo   = Ca->colspi()[0];

    size_t mem = 0.9 * memory_ / 8L;
    // => Target <= //

    auto response = std::make_shared<Matrix>("UHF Response",3*natom,3*natom);

    // => Response Utility File <= //

    psio_->open(PSIF_HESS,PSIO_OPEN_NEW);

    // =>  These functions write alpha/beta components of intermediates to disk <= //
    // Overlap derivatives
    overlap_deriv(Ca, Ca_occ, Ca_vir, nso, naocc, navir, true);
    overlap_deriv(Cb, Cb_occ, Cb_vir, nso, nbocc, nbvir, false);

    // Kinetic derivatives
    kinetic_deriv(Ca, Ca_occ, nso, naocc, navir, true);
    kinetic_deriv(Cb, Cb_occ, nso, nbocc, nbvir, false);

    // Effective core potential derivatives (Epi)
    if (basisset_->has_ECP()) {
#ifdef USING_ecpint
        ecp_deriv(Ca, Ca_occ, nso, naocc, navir, true);
        ecp_deriv(Cb, Cb_occ, nso, nbocc, nbvir, false);
#endif
    }

    // Potential derivatives (Vpi)
    potential_deriv(Ca, Ca_occ, nso, naocc, navir, true);
    potential_deriv(Cb, Cb_occ, nso, nbocc, nbvir, false);

    // Jpi/Kpi
    JK_deriv1(Da, Ca, Ca_occ, Db, nso, naocc, navir, true);
    JK_deriv1(Db, Cb, Cb_occ, Da, nso, nbocc, nbvir, false);

    std::shared_ptr<JK> jk;
    jk = JK::build_JK(basisset_, get_basisset("DF_BASIS_SCF"), options_, false, mem);

    jk->set_memory(mem);
    jk->initialize();

    // Jpi/Kpi
    JK_deriv2(jk,mem, Ca, Ca_occ, Cb, Cb_occ, nso, naocc, nbocc, navir, true);
    JK_deriv2(jk,mem, Cb, Cb_occ, Ca, Ca_occ, nso, nbocc, naocc, nbvir, false);

    VXC_deriv(Ca, Ca_occ, nso, naocc, navir, true);
    VXC_deriv(Cb, Cb_occ, nso, nbocc, nbvir, false);

    assemble_Fock(naocc, navir,true);
    assemble_Fock(nbocc, nbvir,false);

    assemble_B(eps_aocc, naocc, navir,true);
    assemble_B(eps_bocc, nbocc, nbvir,false);

    // => CPHF (Uai) <= //
    // This needs both alpha and beta components simultaneously,
    // so we'll leave it here rather than in a separate function
    {
        uhf_wfn_->set_jk(jk);

        // using naocc here
        size_t per_A = 3L * nso * nso + 1L * naocc * nso;
        size_t max_A = (mem / 2L) / per_A;
        max_A = (max_A > 3 * natom ? 3 * natom : max_A);

        psio_address next_Baia = PSIO_ZERO;
        psio_address next_Uaia = PSIO_ZERO;
        psio_address next_Baib = PSIO_ZERO;
        psio_address next_Uaib = PSIO_ZERO;

        auto Ta = std::make_shared<Matrix>("Ta",navir,naocc);
        auto Tb = std::make_shared<Matrix>("Tb",nbvir,nbocc);
        double** Tap = Ta->pointer();
        double** Tbp = Tb->pointer();

        for (int A = 0; A < 3 * natom; A+=max_A) {
            int nA = max_A;
            if (A + max_A >= 3 * natom) {
                nA = 3 * natom - A;
            }

            std::vector<SharedMatrix> b_vecs;
            // Fill b
            for (int a = 0; a < nA; a++) {
                std::stringstream ss;
                ss << "Perturbation " << a + A;
                
                // Alpha first
                auto Ba = std::make_shared<Matrix>(ss.str(),naocc,navir);
                psio_->read(PSIF_HESS,"Bai^A_a",(char*)Tap[0], static_cast<size_t> (navir) * naocc * sizeof(double),next_Baia,&next_Baia);
                double** Bap = Ba->pointer();
                for (int i = 0; i < naocc; i++) {
                    C_DCOPY(navir,&Tap[0][i],naocc,Bap[i],1);
                }

                // Then beta
                auto Bb = std::make_shared<Matrix>(ss.str(),nbocc,nbvir);
                psio_->read(PSIF_HESS,"Bai^A_b",(char*)Tbp[0], static_cast<size_t> (nbvir) * nbocc * sizeof(double),next_Baib,&next_Baib);
                double** Bbp = Bb->pointer();
                for (int i = 0; i < nbocc; i++) {
                    C_DCOPY(nbvir,&Tbp[0][i],nbocc,Bbp[i],1);
                }
                
                // U matrices from CPHF will accordingly return in alternating alpha/beta order
                b_vecs.push_back(Ba);
                b_vecs.push_back(Bb);
            }

#ifdef USING_BrianQC
    brianCPHFLeftSideFlag = true;
#endif
            auto u_matrices = uhf_wfn_->cphf_solve(b_vecs, options_.get_double("SOLVER_CONVERGENCE"),
                                                   options_.get_int("SOLVER_MAXITER"), print_);

#ifdef USING_BrianQC
    brianCPHFLeftSideFlag = false;
#endif
            // Result in x
            // Write Uas to disk
            for (int a = 0; a < nA; a++) {
                std::stringstream ss;
                ss << "Perturbation " << a + A;
                u_matrices[2*a]->scale(-1);
                double** Xp = u_matrices[2*a]->pointer();
                for (int i = 0; i < naocc; i++) {
                    C_DCOPY(navir,Xp[i],1,&Tap[0][i],naocc);
                }
                psio_->write(PSIF_HESS,"Uai^A_a",(char*)Tap[0], static_cast<size_t> (navir) * naocc * sizeof(double),next_Uaia,&next_Uaia);
            }
            // Write Ubs to disk
            for (int a = 0; a < nA; a++) {
                std::stringstream ss;
                ss << "Perturbation " << a + A;
                u_matrices[2*a+1]->scale(-1);
                double** Xp = u_matrices[2*a+1]->pointer();
                // U matrices stored A/B
                for (int i = 0; i < nbocc; i++) {
                    C_DCOPY(nbvir,Xp[i],1,&Tbp[0][i],nbocc);
                }
                psio_->write(PSIF_HESS,"Uai^A_b",(char*)Tbp[0], static_cast<size_t> (nbvir) * nbocc * sizeof(double),next_Uaib,&next_Uaib);
            }
        }
    }


    assemble_U(naocc, navir, true);
    assemble_U(nbocc, nbvir, false);

    assemble_Q(jk, Ca, Ca_occ, Cb, Cb_occ, nso, naocc, nbocc, navir, true);
    assemble_Q(jk, Cb, Cb_occ, Ca, Ca_occ, nso, nbocc, naocc, nbvir, false);
    jk.reset();

    // => Zipper <= //
    {
        double* eao_p = eps_aocc->pointer();
        double* ebo_p = eps_bocc->pointer();

        double* ea_p = eps_a->pointer();
        double* eb_p = eps_b->pointer();

        size_t memory = 0.9 * memory_ / 8L;
        size_t npia = nmo * (size_t) naocc;
        size_t npib = nmo * (size_t) nbocc;
        size_t max_a = memory / (3L * npia);
        size_t max_b = memory / (3L * npib);
        max_a = (max_a > 3 * natom ? 3 * natom : max_a);
        max_b = (max_b > 3 * natom ? 3 * natom : max_b);

        auto La = std::make_shared<Matrix>("La",max_a * nmo, naocc);
        auto Ra = std::make_shared<Matrix>("Ra",max_a * nmo, naocc);
        auto Ta = std::make_shared<Matrix>("Ta",max_a * nmo, naocc);
        double** Lap = La->pointer();
        double** Rap = Ra->pointer();
        double** Tap = Ta->pointer();

        auto Lb = std::make_shared<Matrix>("Lb",max_b * nmo, nbocc);
        auto Rb = std::make_shared<Matrix>("Rb",max_b * nmo, nbocc);
        auto Tb = std::make_shared<Matrix>("Tb",max_b * nmo, nbocc);
        double** Lbp = Lb->pointer();
        double** Rbp = Rb->pointer();
        double** Tbp = Tb->pointer();

        double** Hp = response->pointer();

        // U^A F^B
        for (int A = 0; A < 3 * natom; A+=max_a) {
            int nA = (A + max_a >= 3 * natom ? 3 * natom - A : max_a);
            psio_address nextA = psio_get_address(PSIO_ZERO, A * npia * sizeof(double));
            psio_->read(PSIF_HESS,"Upi^A_a",(char*)Lap[0],sizeof(double) * nA * npia,nextA,&nextA);
            for (int B = 0; B < 3 * natom; B+=max_a) {
                int nB = (B + max_a >= 3 * natom ? 3 * natom - B : max_a);
                psio_address nextB = psio_get_address(PSIO_ZERO, B * npia * sizeof(double));
                psio_->read(PSIF_HESS,"Fpi^A_a",(char*)Rap[0],sizeof(double) * nB * npia,nextB,&nextB);
                for (int a = 0; a < nA; a++) {
                    for (int b = 0; b < nB; b++) {
                        Hp[A + a][B + b] += 2.0 * C_DDOT(npia,Lap[0] + a * npia,1,Rap[0] + b * npia,1);
                    }
                }
            }
        }
        for (int A = 0; A < 3 * natom; A+=max_b) {
            int nA = (A + max_b >= 3 * natom ? 3 * natom - A : max_b);
            psio_address nextA = psio_get_address(PSIO_ZERO, A * npib * sizeof(double));
            psio_->read(PSIF_HESS,"Upi^A_b",(char*)Lbp[0],sizeof(double) * nA * npib,nextA,&nextA);
            for (int B = 0; B < 3 * natom; B+=max_b) {
                int nB = (B + max_b >= 3 * natom ? 3 * natom - B : max_b);
                psio_address nextB = psio_get_address(PSIO_ZERO, B * npib * sizeof(double));
                psio_->read(PSIF_HESS,"Fpi^A_b",(char*)Rbp[0],sizeof(double) * nB * npib,nextB,&nextB);
                for (int a = 0; a < nA; a++) {
                    for (int b = 0; b < nB; b++) {
                        Hp[A + a][B + b] += 2.0 * C_DDOT(npib,Lbp[0] + a * npib,1,Rbp[0] + b * npib,1);
                    }
                }
            }
        }

        // F^A U^B
        for (int A = 0; A < 3 * natom; A+=max_a) {
            int nA = (A + max_a >= 3 * natom ? 3 * natom - A : max_a);
            psio_address nextA = psio_get_address(PSIO_ZERO, A * npia * sizeof(double));
            psio_->read(PSIF_HESS,"Fpi^A_a",(char*)Lap[0],sizeof(double) * nA * npia,nextA,&nextA);
            for (int B = 0; B < 3 * natom; B+=max_a) {
                int nB = (B + max_a >= 3 * natom ? 3 * natom - B : max_a);
                psio_address nextB = psio_get_address(PSIO_ZERO, B * npia * sizeof(double));
                psio_->read(PSIF_HESS,"Upi^A_a",(char*)Rap[0],sizeof(double) * nB * npia,nextB,&nextB);
                for (int a = 0; a < nA; a++) {
                    for (int b = 0; b < nB; b++) {
                        Hp[A + a][B + b] += 2.0 * C_DDOT(npia,Lap[0] + a * npia,1,Rap[0] + b * npia,1);
                    }
                }
            }
        }
        for (int A = 0; A < 3 * natom; A+=max_b) {
            int nA = (A + max_b >= 3 * natom ? 3 * natom - A : max_b);
            psio_address nextA = psio_get_address(PSIO_ZERO, A * npib * sizeof(double));
            psio_->read(PSIF_HESS,"Fpi^A_b",(char*)Lbp[0],sizeof(double) * nA * npib,nextA,&nextA);
            for (int B = 0; B < 3 * natom; B+=max_b) {
                int nB = (B + max_b >= 3 * natom ? 3 * natom - B : max_b);
                psio_address nextB = psio_get_address(PSIO_ZERO, B * npib * sizeof(double));
                psio_->read(PSIF_HESS,"Upi^A_b",(char*)Rbp[0],sizeof(double) * nB * npib,nextB,&nextB);
                for (int a = 0; a < nA; a++) {
                    for (int b = 0; b < nB; b++) {
                        Hp[A + a][B + b] += 2.0 * C_DDOT(npib,Lbp[0] + a * npib,1,Rbp[0] + b * npib,1);
                    }
                }
            }
        }


        // U^A U^B
        // N.B. We use the relationship U^a_ia = -U^a_ai - S^a_ai
        psio_address junk;
        for (int A = 0; A < 3 * natom; A+=max_a) {
            int nA = (A + max_a >= 3 * natom ? 3 * natom - A : max_a);
            psio_address nextA = psio_get_address(PSIO_ZERO, A * npia * sizeof(double));
            psio_->read(PSIF_HESS,"Upi^A_a",(char*)Lap[0],sizeof(double) * nA * npia,nextA,&junk);
            psio_->read(PSIF_HESS,"Spi^A_a",(char*)Tap[0],sizeof(double) * nA * npia,nextA,&junk);
            La->add(Ta);
            for (int i = 0; i < naocc; i++)
                C_DSCAL(static_cast<size_t> (nA)*nmo,eao_p[i],&Lap[0][i],naocc);
            for (int B = 0; B < 3 * natom; B+=max_a) {
                int nB = (B + max_a >= 3 * natom ? 3 * natom - B : max_a);
                psio_address nextB = psio_get_address(PSIO_ZERO, B * npia * sizeof(double));
                psio_->read(PSIF_HESS,"Upi^A_a",(char*)Rap[0],sizeof(double) * nB * npia,nextB,&junk);
                psio_->read(PSIF_HESS,"Spi^A_a",(char*)Tap[0],sizeof(double) * nB * npia,nextB,&junk);
                Ra->add(Ta);
                for (int a = 0; a < nA; a++) {
                    for (int b = 0; b < nB; b++) {
                        Hp[A + a][B + b] -= C_DDOT(npia,Lap[0] + a * npia,1,Rap[0] + b * npia,1);
                        Hp[B + b][A + a] -= C_DDOT(npia,Lap[0] + a * npia,1,Rap[0] + b * npia,1);
                    }
                }
            }
        }
        for (int A = 0; A < 3 * natom; A+=max_b) {
            int nA = (A + max_b >= 3 * natom ? 3 * natom - A : max_b);
            psio_address nextA = psio_get_address(PSIO_ZERO, A * npib * sizeof(double));
            psio_->read(PSIF_HESS,"Upi^A_b",(char*)Lbp[0],sizeof(double) * nA * npib,nextA,&junk);
            psio_->read(PSIF_HESS,"Spi^A_b",(char*)Tbp[0],sizeof(double) * nA * npib,nextA,&junk);
            Lb->add(Tb);
            for (int i = 0; i < nbocc; i++)
                C_DSCAL(static_cast<size_t> (nA)*nmo,ebo_p[i],&Lbp[0][i],nbocc);
            for (int B = 0; B < 3 * natom; B+=max_b) {
                int nB = (B + max_b >= 3 * natom ? 3 * natom - B : max_b);
                psio_address nextB = psio_get_address(PSIO_ZERO, B * npib * sizeof(double));
                psio_->read(PSIF_HESS,"Upi^A_b",(char*)Rbp[0],sizeof(double) * nB * npib,nextB,&junk);
                psio_->read(PSIF_HESS,"Spi^A_b",(char*)Tbp[0],sizeof(double) * nB * npib,nextB,&junk);
                Rb->add(Tb);
                for (int a = 0; a < nA; a++) {
                    for (int b = 0; b < nB; b++) {
                        Hp[A + a][B + b] -= C_DDOT(npib,Lbp[0] + a * npib,1,Rbp[0] + b * npib,1);
                        Hp[B + b][A + a] -= C_DDOT(npib,Lbp[0] + a * npib,1,Rbp[0] + b * npib,1);
                    }
                }
            }
        }
        // S^A S^B
        for (int A = 0; A < 3 * natom; A+=max_a) {
            int nA = (A + max_a >= 3 * natom ? 3 * natom - A : max_a);
            psio_address nextA = psio_get_address(PSIO_ZERO, A * npia * sizeof(double));
            psio_->read(PSIF_HESS,"Spi^A_a",(char*)Lap[0],sizeof(double) * nA * npia,nextA,&nextA);
            for (int i = 0; i < naocc; i++)
                C_DSCAL(static_cast<size_t> (nA)*nmo,eao_p[i],&Lap[0][i],naocc);
            for (int B = 0; B < 3 * natom; B+=max_a) {
                int nB = (B + max_a >= 3 * natom ? 3 * natom - B : max_a);
                psio_address nextB = psio_get_address(PSIO_ZERO, B * npia * sizeof(double));
                psio_->read(PSIF_HESS,"Spi^A_a",(char*)Rap[0],sizeof(double) * nB * npia,nextB,&nextB);
                for (int a = 0; a < nA; a++) {
                    for (int b = 0; b < nB; b++) {
                        Hp[A + a][B + b] += C_DDOT(npia,Lap[0] + a * npia,1,Rap[0] + b * npia,1);
                        Hp[B + b][A + a] += C_DDOT(npia,Lap[0] + a * npia,1,Rap[0] + b * npia,1);
                    }
                }
            }
        }
        for (int A = 0; A < 3 * natom; A+=max_b) {
            int nA = (A + max_b >= 3 * natom ? 3 * natom - A : max_b);
            psio_address nextA = psio_get_address(PSIO_ZERO, A * npib * sizeof(double));
            psio_->read(PSIF_HESS,"Spi^A_b",(char*)Lbp[0],sizeof(double) * nA * npib,nextA,&nextA);
            for (int i = 0; i < nbocc; i++)
                C_DSCAL(static_cast<size_t> (nA)*nmo,ebo_p[i],&Lbp[0][i],nbocc);
            for (int B = 0; B < 3 * natom; B+=max_b) {
                int nB = (B + max_b >= 3 * natom ? 3 * natom - B : max_b);
                psio_address nextB = psio_get_address(PSIO_ZERO, B * npib * sizeof(double));
                psio_->read(PSIF_HESS,"Spi^A_b",(char*)Rbp[0],sizeof(double) * nB * npib,nextB,&nextB);
                for (int a = 0; a < nA; a++) {
                    for (int b = 0; b < nB; b++) {
                        Hp[A + a][B + b] += C_DDOT(npib,Lbp[0] + a * npib,1,Rbp[0] + b * npib,1);
                        Hp[B + b][A + a] += C_DDOT(npib,Lbp[0] + a * npib,1,Rbp[0] + b * npib,1);
                    }
                }
            }
        }

        // U^A U^B \epsilon
        for (int A = 0; A < 3 * natom; A+=max_a) {
            int nA = (A + max_a >= 3 * natom ? 3 * natom - A : max_a);
            psio_address nextA = psio_get_address(PSIO_ZERO, A * npia * sizeof(double));
            psio_->read(PSIF_HESS,"Upi^A_a",(char*)Lap[0],sizeof(double) * nA * npia,nextA,&nextA);
            double* Tap = Lap[0];
            for (int a = 0; a < nA; a++) {
                for (int p = 0; p < nmo; p++) {
                    C_DSCAL(naocc,ea_p[p],Tap,1);
                    Tap += naocc;
                }
            }
            for (int B = 0; B < 3 * natom; B+=max_a) {
                int nB = (B + max_a >= 3 * natom ? 3 * natom - B : max_a);
                psio_address nextB = psio_get_address(PSIO_ZERO, B * npia * sizeof(double));
                psio_->read(PSIF_HESS,"Upi^A_a",(char*)Rap[0],sizeof(double) * nB * npia,nextB,&nextB);
                for (int a = 0; a < nA; a++) {
                    for (int b = 0; b < nB; b++) {
                        Hp[A + a][B + b] += 2.0 * C_DDOT(npia,Lap[0] + a * npia,1,Rap[0] + b * npia,1);
                    }
                }
            }
        }
        for (int A = 0; A < 3 * natom; A+=max_b) {
            int nA = (A + max_b >= 3 * natom ? 3 * natom - A : max_b);
            psio_address nextA = psio_get_address(PSIO_ZERO, A * npib * sizeof(double));
            psio_->read(PSIF_HESS,"Upi^A_b",(char*)Lbp[0],sizeof(double) * nA * npib,nextA,&nextA);
            double* Tbp = Lbp[0];
            for (int a = 0; a < nA; a++) {
                for (int p = 0; p < nmo; p++) {
                    C_DSCAL(nbocc,eb_p[p],Tbp,1);
                    Tbp += nbocc;
                }
            }
            for (int B = 0; B < 3 * natom; B+=max_b) {
                int nB = (B + max_b >= 3 * natom ? 3 * natom - B : max_b);
                psio_address nextB = psio_get_address(PSIO_ZERO, B * npib * sizeof(double));
                psio_->read(PSIF_HESS,"Upi^A_b",(char*)Rbp[0],sizeof(double) * nB * npib,nextB,&nextB);
                for (int a = 0; a < nA; a++) {
                    for (int b = 0; b < nB; b++) {
                        Hp[A + a][B + b] += 2.0 * C_DDOT(npib,Lbp[0] + a * npib,1,Rbp[0] + b * npib,1);
                    }
                }
            }
        }

        // U^A Q^B
        for (int A = 0; A < 3 * natom; A+=max_a) {
            int nA = (A + max_a >= 3 * natom ? 3 * natom - A : max_a);
            psio_address nextA = psio_get_address(PSIO_ZERO, A * npia * sizeof(double));
            psio_->read(PSIF_HESS,"Upi^A_a",(char*)Lap[0],sizeof(double) * nA * npia,nextA,&nextA);
            for (int B = 0; B < 3 * natom; B+=max_a) {
                int nB = (B + max_a >= 3 * natom ? 3 * natom - B : max_a);
                psio_address nextB = psio_get_address(PSIO_ZERO, B * npia * sizeof(double));
                psio_->read(PSIF_HESS,"Qpi^A_a",(char*)Rap[0],sizeof(double) * nB * npia,nextB,&nextB);
                for (int a = 0; a < nA; a++) {
                    for (int b = 0; b < nB; b++) {
                        Hp[A + a][B + b] += 2.0 * C_DDOT(npia,Lap[0] + a * npia,1,Rap[0] + b * npia,1);
                    }
                }
            }
        }
        for (int A = 0; A < 3 * natom; A+=max_b) {
            int nA = (A + max_b >= 3 * natom ? 3 * natom - A : max_b);
            psio_address nextA = psio_get_address(PSIO_ZERO, A * npib * sizeof(double));
            psio_->read(PSIF_HESS,"Upi^A_b",(char*)Lbp[0],sizeof(double) * nA * npib,nextA,&nextA);
            for (int B = 0; B < 3 * natom; B+=max_b) {
                int nB = (B + max_b >= 3 * natom ? 3 * natom - B : max_b);
                psio_address nextB = psio_get_address(PSIO_ZERO, B * npib * sizeof(double));
                psio_->read(PSIF_HESS,"Qpi^A_b",(char*)Rbp[0],sizeof(double) * nB * npib,nextB,&nextB);
                for (int a = 0; a < nA; a++) {
                    for (int b = 0; b < nB; b++) {
                        Hp[A + a][B + b] += 2.0 * C_DDOT(npib,Lbp[0] + a * npib,1,Rbp[0] + b * npib,1);
                    }
                }
            }
        }


        // Full symmetrization
        for (int A = 0; A < 3 * natom; A++) {
            for (int B = 0; B < 3 * natom; B++) {
                Hp[A][B] = Hp[B][A] = 0.5*(Hp[A][B] + Hp[B][A]);
            }
        }
    }
    psio_->close(PSIF_HESS,0);
    return response;
}


void USCFDeriv::overlap_deriv(std::shared_ptr<Matrix> C, 
                              std::shared_ptr<Matrix> Cocc,
                              std::shared_ptr<Matrix> Cvir,
                              int nso, int nocc, int nvir, bool alpha)
{
    // Overlap derivatives

    double** Cp  = C->pointer();  
    double** Cop = Cocc->pointer();
    double** Cvp = Cvir->pointer(); 

    std::shared_ptr<OneBodyAOInt> Sint(integral_->ao_overlap(1));
    size_t nmo = nocc + nvir;
    int natom = molecule_->natom();

    auto Smix = std::make_shared<Matrix>("Smix",nso,nocc);
    auto Smiy = std::make_shared<Matrix>("Smiy",nso,nocc);
    auto Smiz = std::make_shared<Matrix>("Smiz",nso,nocc);
    double** Smixp = Smix->pointer();
    double** Smiyp = Smiy->pointer();
    double** Smizp = Smiz->pointer();

    auto Sai = std::make_shared<Matrix>("Sai",nvir,nocc);
    double** Saip = Sai->pointer();
    auto Sij = std::make_shared<Matrix>("Sij",nocc,nocc);
    double** Sijp = Sij->pointer();
    auto Spi = std::make_shared<Matrix>("Spi",nmo,nocc);
    double** Spip = Spi->pointer();

    psio_address next_Sai = PSIO_ZERO;
    psio_address next_Sij = PSIO_ZERO;
    psio_address next_Smi = PSIO_ZERO;
    psio_address next_Spi = PSIO_ZERO;
    
#ifdef USING_BrianQC
    brianInt maxSegmentSize;
    brianInt maxSegmentAtomCount;
    brianInt segmentAtomCount;
    brianInt segmentAtomIndexStart;
    
    std::vector<std::shared_ptr<Matrix>> Smnx;
    std::vector<std::shared_ptr<Matrix>> Smny;
    std::vector<std::shared_ptr<Matrix>> Smnz;
    std::vector<double*> Smn;
    
    if (brianEnable) {
        brianCPHFMaxSegmentSize(&brianCookie, &maxSegmentSize);
        maxSegmentAtomCount = maxSegmentSize / 3;
        segmentAtomCount = -1;
        segmentAtomIndexStart = -1;
        
        Smnx.resize(maxSegmentAtomCount);
        Smny.resize(maxSegmentAtomCount);
        Smnz.resize(maxSegmentAtomCount);
        Smn.resize(maxSegmentAtomCount * 3);
        for (int i = 0; i < maxSegmentAtomCount; i++) {
            Smnx[i] = std::make_shared<Matrix>("Smnx", nso, nso);
            Smny[i] = std::make_shared<Matrix>("Smny", nso, nso);
            Smnz[i] = std::make_shared<Matrix>("Smnz", nso, nso);
            Smn[i * 3 + 0] = Smnx[i]->get_pointer();
            Smn[i * 3 + 1] = Smny[i]->get_pointer();
            Smn[i * 3 + 2] = Smnz[i]->get_pointer();
        }
    }
#endif

    // Get the filenames right
    auto Smi_str = (alpha) ? "Smi^A_a" : "Smi^A_b";
    auto Sai_str = (alpha) ? "Sai^A_a" : "Sai^A_b";
    auto Sij_str = (alpha) ? "Sij^A_a" : "Sij^A_b";
    auto Spi_str = (alpha) ? "Spi^A_a" : "Spi^A_b";

    for (int A = 0; A < 3*natom; A++) {
        psio_->write(PSIF_HESS,Smi_str,(char*)Smixp[0], static_cast<size_t> (nso) * nocc * sizeof(double),next_Smi,&next_Smi);
    }
    for (int A = 0; A < 3*natom; A++) {
        psio_->write(PSIF_HESS,Sai_str,(char*)Saip[0], static_cast<size_t> (nvir) * nocc * sizeof(double),next_Sai,&next_Sai);
    }
    for (int A = 0; A < 3*natom; A++) {
        psio_->write(PSIF_HESS,Sij_str,(char*)Sijp[0], static_cast<size_t> (nocc) * nocc * sizeof(double),next_Sij,&next_Sij);
    }
    for (int A = 0; A < 3*natom; A++) {
        psio_->write(PSIF_HESS,Spi_str,(char*)Spip[0], static_cast<size_t> (nmo) * nocc * sizeof(double),next_Spi,&next_Spi);
    }

    next_Smi = PSIO_ZERO;
    next_Sai = PSIO_ZERO;
    next_Sij = PSIO_ZERO;
    next_Spi = PSIO_ZERO;

    for (int A = 0; A < natom; A++) {
#ifdef USING_BrianQC
        if (brianEnable) {
            if (segmentAtomCount < 0 || A < segmentAtomIndexStart || A >= (segmentAtomIndexStart + segmentAtomCount)) {
                segmentAtomIndexStart = A;
                segmentAtomCount = (segmentAtomIndexStart + maxSegmentAtomCount > natom) ? (natom - segmentAtomIndexStart) : maxSegmentAtomCount;
                
                brianInt integralType = BRIAN_INTEGRAL_TYPE_OVERLAP;
                brianCPHFBuild1eDeriv(&brianCookie, &integralType, &segmentAtomCount, &segmentAtomIndexStart, Smn.data());
            }
            
            C_DGEMM('N', 'N', nso, nocc, nso, 2.0, Smnx[A - segmentAtomIndexStart]->get_pointer(), nso, Cocc->get_pointer(), nocc, 0.0, Smix->get_pointer(), nocc);
            C_DGEMM('N', 'N', nso, nocc, nso, 2.0, Smny[A - segmentAtomIndexStart]->get_pointer(), nso, Cocc->get_pointer(), nocc, 0.0, Smiy->get_pointer(), nocc);
            C_DGEMM('N', 'N', nso, nocc, nso, 2.0, Smnz[A - segmentAtomIndexStart]->get_pointer(), nso, Cocc->get_pointer(), nocc, 0.0, Smiz->get_pointer(), nocc);
        } else {
#endif
        Smix->zero();
        Smiy->zero();
        Smiz->zero();
        const auto& shell_pairs = Sint->shellpairs();
        size_t n_pairs = shell_pairs.size();

        for (size_t p = 0; p < n_pairs; ++p) {
            auto P = shell_pairs[p].first;
            auto Q = shell_pairs[p].second;
            const auto &shellP = basisset_->shell(P);
            const auto &shellQ = basisset_->shell(Q);
            int aP = shellP.ncenter();
            int aQ = shellQ.ncenter();
            if ((aP != A && aQ != A) || aP == aQ) continue;
            Sint->compute_shell_deriv1(P, Q);
            const auto &buffers = Sint->buffers();
            int nP = shellP.nfunction();
            int nQ = shellQ.nfunction();
            int oP = shellP.function_index();
            int oQ = shellQ.function_index();
            const double* buffer2;

            double scale = P == Q ? 1.0 : 2.0;
            if (aP == A) {
                // Px
                buffer2 = buffers[0];
                for (int p = 0; p < nP; p++) {
                    for (int q = 0; q < nQ; q++) {
                        C_DAXPY(nocc,scale*(*buffer2),Cop[q + oQ],1,Smixp[p + oP],1);
                        C_DAXPY(nocc,scale*(*buffer2++), Cop[p + oP],1,Smixp[q + oQ],1);
                    }
                }
                // Py
                buffer2 = buffers[1];
                for (int p = 0; p < nP; p++) {
                    for (int q = 0; q < nQ; q++) {
                        C_DAXPY(nocc,scale*(*buffer2),Cop[q + oQ],1,Smiyp[p + oP],1);
                        C_DAXPY(nocc,scale*(*buffer2++), Cop[p + oP],1,Smiyp[q + oQ],1);
                    }
                }
                // Pz
                buffer2 = buffers[2];
                for (int p = 0; p < nP; p++) {
                    for (int q = 0; q < nQ; q++) {
                        C_DAXPY(nocc,scale*(*buffer2),Cop[q + oQ],1,Smizp[p + oP],1);
                        C_DAXPY(nocc,scale*(*buffer2++), Cop[p + oP],1,Smizp[q + oQ],1);
                    }
                }
            }
            if (aQ == A) {
                // Qx
                buffer2 = buffers[3];
                for (int p = 0; p < nP; p++) {
                    for (int q = 0; q < nQ; q++) {
                        C_DAXPY(nocc,scale*(*buffer2),Cop[q + oQ],1,Smixp[p + oP],1);
                        C_DAXPY(nocc,scale*(*buffer2++), Cop[p + oP],1,Smixp[q + oQ],1);
                    }
                }
                // Qy
                buffer2 = buffers[4];
                for (int p = 0; p < nP; p++) {
                    for (int q = 0; q < nQ; q++) {
                        C_DAXPY(nocc,scale*(*buffer2),Cop[q + oQ],1,Smiyp[p + oP],1);
                        C_DAXPY(nocc,scale*(*buffer2++), Cop[p + oP],1,Smiyp[q + oQ],1);
                    }
                }
                // Qz
                buffer2 = buffers[5];
                for (int p = 0; p < nP; p++) {
                    for (int q = 0; q < nQ; q++) {
                        C_DAXPY(nocc,scale*(*buffer2),Cop[q + oQ],1,Smizp[p + oP],1);
                        C_DAXPY(nocc,scale*(*buffer2++), Cop[p + oP],1,Smizp[q + oQ],1);
                    }
                }
            }
        }
#ifdef USING_BrianQC
        }
#endif
        // Smi_x
        psio_->write(PSIF_HESS,Smi_str,(char*)Smixp[0], static_cast<size_t> (nso) * nocc * sizeof(double),next_Smi,&next_Smi);
        // Smi_y
        psio_->write(PSIF_HESS,Smi_str,(char*)Smiyp[0], static_cast<size_t> (nso) * nocc * sizeof(double),next_Smi,&next_Smi);
        // Smi_z
        psio_->write(PSIF_HESS,Smi_str,(char*)Smizp[0], static_cast<size_t> (nso) * nocc * sizeof(double),next_Smi,&next_Smi);

        // Sai_x
        C_DGEMM('T','N',nvir,nocc,nso,0.5,Cvp[0],nvir,Smixp[0],nocc,0.0,Saip[0],nocc);
        psio_->write(PSIF_HESS,Sai_str,(char*)Saip[0], static_cast<size_t> (nvir) * nocc * sizeof(double),next_Sai,&next_Sai);
        // Sai_y
        C_DGEMM('T','N',nvir,nocc,nso,0.5,Cvp[0],nvir,Smiyp[0],nocc,0.0,Saip[0],nocc);
        psio_->write(PSIF_HESS,Sai_str,(char*)Saip[0], static_cast<size_t> (nvir) * nocc * sizeof(double),next_Sai,&next_Sai);
        // Sai_z
        C_DGEMM('T','N',nvir,nocc,nso,0.5,Cvp[0],nvir,Smizp[0],nocc,0.0,Saip[0],nocc);
        psio_->write(PSIF_HESS,Sai_str,(char*)Saip[0], static_cast<size_t> (nvir) * nocc * sizeof(double),next_Sai,&next_Sai);

        // Sij_x
        C_DGEMM('T','N',nocc,nocc,nso,0.5,Cop[0],nocc,Smixp[0],nocc,0.0,Sijp[0],nocc);
        psio_->write(PSIF_HESS,Sij_str,(char*)Sijp[0], static_cast<size_t> (nocc) * nocc * sizeof(double),next_Sij,&next_Sij);
        // Sij_y
        C_DGEMM('T','N',nocc,nocc,nso,0.5,Cop[0],nocc,Smiyp[0],nocc,0.0,Sijp[0],nocc);
        psio_->write(PSIF_HESS,Sij_str,(char*)Sijp[0], static_cast<size_t> (nocc) * nocc * sizeof(double),next_Sij,&next_Sij);
        // Sij_z
        C_DGEMM('T','N',nocc,nocc,nso,0.5,Cop[0],nocc,Smizp[0],nocc,0.0,Sijp[0],nocc);
        psio_->write(PSIF_HESS,Sij_str,(char*)Sijp[0], static_cast<size_t> (nocc) * nocc * sizeof(double),next_Sij,&next_Sij);

        // Spi_x
        C_DGEMM('T','N',nmo,nocc,nso,0.5,Cp[0],nmo,Smixp[0],nocc,0.0,Spip[0],nocc);
        psio_->write(PSIF_HESS,Spi_str,(char*)Spip[0], static_cast<size_t> (nmo) * nocc * sizeof(double),next_Spi,&next_Spi);
        // Spi_y
        C_DGEMM('T','N',nmo,nocc,nso,0.5,Cp[0],nmo,Smiyp[0],nocc,0.0,Spip[0],nocc);
        psio_->write(PSIF_HESS,Spi_str,(char*)Spip[0], static_cast<size_t> (nmo) * nocc * sizeof(double),next_Spi,&next_Spi);
        // Spi_z
        C_DGEMM('T','N',nmo,nocc,nso,0.5,Cp[0],nmo,Smizp[0],nocc,0.0,Spip[0],nocc);
        psio_->write(PSIF_HESS,Spi_str,(char*)Spip[0], static_cast<size_t> (nmo) * nocc * sizeof(double),next_Spi,&next_Spi);
    }
}

void USCFDeriv::kinetic_deriv(std::shared_ptr<Matrix> C, 
                              std::shared_ptr<Matrix> Cocc,
                              int nso, int nocc, int nvir, bool alpha)
{
    // Kinetic derivatives
    std::shared_ptr<OneBodyAOInt> Tint(integral_->ao_kinetic(1));
    size_t nmo = nocc + nvir;
    int natom = molecule_->natom();

    double** Cp  = C->pointer();  
    double** Cop = Cocc->pointer();

    auto Tmix = std::make_shared<Matrix>("Tmix",nso,nocc);
    auto Tmiy = std::make_shared<Matrix>("Tmiy",nso,nocc);
    auto Tmiz = std::make_shared<Matrix>("Tmiz",nso,nocc);
    double** Tmixp = Tmix->pointer();
    double** Tmiyp = Tmiy->pointer();
    double** Tmizp = Tmiz->pointer();

    auto Tpi = std::make_shared<Matrix>("Tpi",nmo,nocc);
    double** Tpip = Tpi->pointer();
    psio_address next_Tpi = PSIO_ZERO;

#ifdef USING_BrianQC
    brianInt maxSegmentSize;
    brianInt maxSegmentAtomCount;
    brianInt segmentAtomCount;
    brianInt segmentAtomIndexStart;
    
    std::vector<std::shared_ptr<Matrix>> Tmnx;
    std::vector<std::shared_ptr<Matrix>> Tmny;
    std::vector<std::shared_ptr<Matrix>> Tmnz;
    std::vector<double*> Tmn;
    
    if (brianEnable) {
        brianCPHFMaxSegmentSize(&brianCookie, &maxSegmentSize);
        maxSegmentAtomCount = maxSegmentSize / 3;
        segmentAtomCount = -1;
        segmentAtomIndexStart = -1;
        
        Tmnx.resize(maxSegmentAtomCount);
        Tmny.resize(maxSegmentAtomCount);
        Tmnz.resize(maxSegmentAtomCount);
        Tmn.resize(maxSegmentAtomCount * 3);
        for (int i = 0; i < maxSegmentAtomCount; i++) {
            Tmnx[i] = std::make_shared<Matrix>("Tmnx", nso, nso);
            Tmny[i] = std::make_shared<Matrix>("Tmny", nso, nso);
            Tmnz[i] = std::make_shared<Matrix>("Tmnz", nso, nso);
            Tmn[i * 3 + 0] = Tmnx[i]->get_pointer();
            Tmn[i * 3 + 1] = Tmny[i]->get_pointer();
            Tmn[i * 3 + 2] = Tmnz[i]->get_pointer();
        }
    }
#endif

    for (int A = 0; A < natom; A++) {
#ifdef USING_BrianQC
        if (brianEnable) {
            if (segmentAtomCount < 0 || A < segmentAtomIndexStart || A >= (segmentAtomIndexStart + segmentAtomCount)) {
                segmentAtomIndexStart = A;
                segmentAtomCount = (segmentAtomIndexStart + maxSegmentAtomCount > natom) ? (natom - segmentAtomIndexStart) : maxSegmentAtomCount;
                
                brianInt integralType = BRIAN_INTEGRAL_TYPE_KINETIC;
                brianCPHFBuild1eDeriv(&brianCookie, &integralType, &segmentAtomCount, &segmentAtomIndexStart, Tmn.data());
            }
            
            C_DGEMM('N', 'N', nso, nocc, nso, 2.0, Tmnx[A - segmentAtomIndexStart]->get_pointer(), nso, Cocc->get_pointer(), nocc, 0.0, Tmix->get_pointer(), nocc);
            C_DGEMM('N', 'N', nso, nocc, nso, 2.0, Tmny[A - segmentAtomIndexStart]->get_pointer(), nso, Cocc->get_pointer(), nocc, 0.0, Tmiy->get_pointer(), nocc);
            C_DGEMM('N', 'N', nso, nocc, nso, 2.0, Tmnz[A - segmentAtomIndexStart]->get_pointer(), nso, Cocc->get_pointer(), nocc, 0.0, Tmiz->get_pointer(), nocc);
        } else {
#endif
        Tmix->zero();
        Tmiy->zero();
        Tmiz->zero();
        const auto& shell_pairs = Tint->shellpairs();
        size_t n_pairs = shell_pairs.size();

        for (size_t p = 0; p < n_pairs; ++p) {
            auto P = shell_pairs[p].first;
            auto Q = shell_pairs[p].second;
            const auto & shellP = basisset_->shell(P);
            const auto & shellQ = basisset_->shell(Q);
            int aP = shellP.ncenter();
            int aQ = shellQ.ncenter();
            if ((aP != A && aQ != A) || aP == aQ) continue;
            Tint->compute_shell_deriv1(P, Q);
            const auto &buffers = Tint->buffers();
            int nP = shellP.nfunction();
            int nQ = shellQ.nfunction();
            int oP = shellP.function_index();
            int oQ = shellQ.function_index();
            const double* buffer2;
            double scale = P == Q ? 1.0 : 2.0;
            if (aP == A) {
                // Px
                buffer2 = buffers[0];
                for (int p = 0; p < nP; p++) {
                    for (int q = 0; q < nQ; q++) {
                        C_DAXPY(nocc,scale*(*buffer2),Cop[q + oQ],1,Tmixp[p + oP],1);
                        C_DAXPY(nocc,scale*(*buffer2++),Cop[p + oP],1,Tmixp[q + oQ],1);
                    }
                }
                // Py
                buffer2 = buffers[1];
                for (int p = 0; p < nP; p++) {
                    for (int q = 0; q < nQ; q++) {
                        C_DAXPY(nocc,scale*(*buffer2),Cop[q + oQ],1,Tmiyp[p + oP],1);
                        C_DAXPY(nocc,scale*(*buffer2++),Cop[p + oP],1,Tmiyp[q + oQ],1);
                    }
                }
                // Pz
                buffer2 = buffers[2];
                for (int p = 0; p < nP; p++) {
                    for (int q = 0; q < nQ; q++) {
                        C_DAXPY(nocc,scale*(*buffer2),Cop[q + oQ],1,Tmizp[p + oP],1);
                        C_DAXPY(nocc,scale*(*buffer2++),Cop[p + oP],1,Tmizp[q + oQ],1);
                    }
                }
            }
            if (aQ == A) {
                // Qx
                buffer2 = buffers[3];
                for (int p = 0; p < nP; p++) {
                    for (int q = 0; q < nQ; q++) {
                        C_DAXPY(nocc,scale*(*buffer2),Cop[q + oQ],1,Tmixp[p + oP],1);
                        C_DAXPY(nocc,scale*(*buffer2++),Cop[p + oP],1,Tmixp[q + oQ],1);
                    }
                }
                // Qy
                buffer2 = buffers[4];
                for (int p = 0; p < nP; p++) {
                    for (int q = 0; q < nQ; q++) {
                        C_DAXPY(nocc,scale*(*buffer2),Cop[q + oQ],1,Tmiyp[p + oP],1);
                        C_DAXPY(nocc,scale*(*buffer2++),Cop[p + oP],1,Tmiyp[q + oQ],1);
                    }
                }
                // Qz
                buffer2 = buffers[5];
                for (int p = 0; p < nP; p++) {
                    for (int q = 0; q < nQ; q++) {
                        C_DAXPY(nocc,scale*(*buffer2),Cop[q + oQ],1,Tmizp[p + oP],1);
                        C_DAXPY(nocc,scale*(*buffer2++),Cop[p + oP],1,Tmizp[q + oQ],1);
                    }
                }
            }
        }
#ifdef USING_BrianQC
        }
#endif

        auto Tpi_str = (alpha) ? "Tpi^A_a" : "Tpi^A_b";

        // Tpi_x
        C_DGEMM('T','N',nmo,nocc,nso,0.5,Cp[0],nmo,Tmixp[0],nocc,0.0,Tpip[0],nocc);
        psio_->write(PSIF_HESS,Tpi_str,(char*)Tpip[0], static_cast<size_t> (nmo) * nocc * sizeof(double),next_Tpi,&next_Tpi);
        // Tpi_y
        C_DGEMM('T','N',nmo,nocc,nso,0.5,Cp[0],nmo,Tmiyp[0],nocc,0.0,Tpip[0],nocc);
        psio_->write(PSIF_HESS,Tpi_str,(char*)Tpip[0], static_cast<size_t> (nmo) * nocc * sizeof(double),next_Tpi,&next_Tpi);
        // Tpi_z
        C_DGEMM('T','N',nmo,nocc,nso,0.5,Cp[0],nmo,Tmizp[0],nocc,0.0,Tpip[0],nocc);
        psio_->write(PSIF_HESS,Tpi_str,(char*)Tpip[0], static_cast<size_t> (nmo) * nocc * sizeof(double),next_Tpi,&next_Tpi);
    }
}

#ifdef USING_ecpint
void USCFDeriv::ecp_deriv(std::shared_ptr<Matrix> C, 
                          std::shared_ptr<Matrix> Cocc,
                          int nso, int nocc, int nvir, bool alpha)
{
    std::shared_ptr<ECPInt> ecpint(dynamic_cast<ECPInt*>(integral_->ao_ecp(1).release()));
    size_t nmo = nocc + nvir;
    int natom = molecule_->natom();

    double** Cp  = C->pointer();  
    double** Cop = Cocc->pointer();

    auto Emix = std::make_shared<Matrix>("Emix", nso, nocc);
    auto Emiy = std::make_shared<Matrix>("Emiy", nso, nocc);
    auto Emiz = std::make_shared<Matrix>("Emiz", nso, nocc);
    double** Emixp = Emix->pointer();
    double** Emiyp = Emiy->pointer();
    double** Emizp = Emiz->pointer();

    // Make a list of all ECP centers
    std::set<int> ecp_centers;
    for (int ecp_shell = 0; ecp_shell < basisset_->n_ecp_shell(); ++ecp_shell){
        const GaussianShell &ecp = basisset_->ecp_shell(ecp_shell);
        ecp_centers.insert(ecp.ncenter());
    }

    auto Epi = std::make_shared<Matrix>("Epi", nmo, nocc);
    double** Epip = Epi->pointer();
    psio_address next_Epi = PSIO_ZERO;


    for (int A = 0; A < natom; A++) {
        Emix->zero();
        Emiy->zero();
        Emiz->zero();
        const auto& shell_pairs = ecpint->shellpairs();
        size_t n_pairs = shell_pairs.size();

        for (size_t p = 0; p < n_pairs; ++p) {
            auto P = shell_pairs[p].first;
            auto Q = shell_pairs[p].second;
            const auto & shellP = basisset_->shell(P);
            const auto & shellQ = basisset_->shell(Q);
            int aP = shellP.ncenter();
            int aQ = shellQ.ncenter();

            // Make a list of all ECP centers and the current basis function pair
            std::set<int> all_centers(ecp_centers.begin(), ecp_centers.end());
            all_centers.insert(aP);
            all_centers.insert(aQ);
            if (all_centers.find(A) == all_centers.end()) continue;

            ecpint->compute_shell_deriv1(P, Q);
            const auto &buffers = ecpint->buffers();
            int nP = shellP.nfunction();
            int nQ = shellQ.nfunction();
            int oP = shellP.function_index();
            int oQ = shellQ.function_index();
            const double* buffer2;

            for (const int center : all_centers) {
                double scale = P == Q ? 1.0 : 2.0;
                if (center == A) {
                    // x
                    buffer2 = buffers[3*center+0];
                    for (int p = 0; p < nP; p++) {
                        for (int q = 0; q < nQ; q++) {
                            C_DAXPY(nocc, scale * (*buffer2), Cop[q + oQ], 1, Emixp[p + oP], 1);
                            C_DAXPY(nocc, scale * (*buffer2++), Cop[p + oP], 1, Emixp[q + oQ], 1);
                        }
                    }
                    // y
                    buffer2 = buffers[3*center+1];
                    for (int p = 0; p < nP; p++) {
                        for (int q = 0; q < nQ; q++) {
                            C_DAXPY(nocc, scale * (*buffer2), Cop[q + oQ], 1, Emiyp[p + oP], 1);
                            C_DAXPY(nocc, scale * (*buffer2++), Cop[p + oP], 1, Emiyp[q + oQ], 1);
                        }
                    }
                    // z
                    buffer2 = buffers[3*center+2];
                    for (int p = 0; p < nP; p++) {
                        for (int q = 0; q < nQ; q++) {
                            C_DAXPY(nocc, scale * (*buffer2), Cop[q + oQ], 1, Emizp[p + oP], 1);
                            C_DAXPY(nocc, scale * (*buffer2++), Cop[p + oP], 1, Emizp[q + oQ], 1);
                        }
                    }
                }
            }
        }
        auto Epi_str = (alpha) ? "Epi^A_a" : "Epi^A_b";
        // Epi_x
        C_DGEMM('T', 'N', nmo, nocc, nso, 0.5, Cp[0], nmo, Emixp[0], nocc, 0.0, Epip[0], nocc);
        psio_->write(PSIF_HESS, Epi_str, (char*)Epip[0], static_cast<size_t>(nmo) * nocc * sizeof(double), next_Epi,
                     &next_Epi);
        // Epi_y
        C_DGEMM('T', 'N', nmo, nocc, nso, 0.5, Cp[0], nmo, Emiyp[0], nocc, 0.0, Epip[0], nocc);
        psio_->write(PSIF_HESS, Epi_str, (char*)Epip[0], static_cast<size_t>(nmo) * nocc * sizeof(double), next_Epi,
                     &next_Epi);
        // Epi_z
        C_DGEMM('T', 'N', nmo, nocc, nso, 0.5, Cp[0], nmo, Emizp[0], nocc, 0.0, Epip[0], nocc);
        psio_->write(PSIF_HESS, Epi_str, (char*)Epip[0], static_cast<size_t>(nmo) * nocc * sizeof(double), next_Epi,
                     &next_Epi);
    }
}
#endif

void USCFDeriv::potential_deriv(std::shared_ptr<Matrix> C, 
                                std::shared_ptr<Matrix> Cocc,
                                int nso, int nocc, int nvir, bool alpha)
{
    std::shared_ptr<OneBodyAOInt> Vint(integral_->ao_potential(1));
    size_t nmo = nocc + nvir;
    int natom = molecule_->natom();

    double** Cp  = C->pointer();  
    double** Cop = Cocc->pointer();

    auto Vmix = std::make_shared<Matrix>("Vmix",nso,nocc);
    auto Vmiy = std::make_shared<Matrix>("Vmiy",nso,nocc);
    auto Vmiz = std::make_shared<Matrix>("Vmiz",nso,nocc);
    double** Vmixp = Vmix->pointer();
    double** Vmiyp = Vmiy->pointer();
    double** Vmizp = Vmiz->pointer();

    auto Vpi = std::make_shared<Matrix>("Vpi",nmo,nocc);
    double** Vpip = Vpi->pointer();
    psio_address next_Vpi = PSIO_ZERO;

#ifdef USING_BrianQC
    brianInt maxSegmentSize;
    brianInt maxSegmentAtomCount;
    brianInt segmentAtomCount;
    brianInt segmentAtomIndexStart;
    
    std::vector<std::shared_ptr<Matrix>> Vmnx;
    std::vector<std::shared_ptr<Matrix>> Vmny;
    std::vector<std::shared_ptr<Matrix>> Vmnz;
    std::vector<double*> Vmn;
    
    if (brianEnable) {
        brianCPHFMaxSegmentSize(&brianCookie, &maxSegmentSize);
        maxSegmentAtomCount = maxSegmentSize / 3;
        segmentAtomCount = -1;
        segmentAtomIndexStart = -1;
        
        Vmnx.resize(maxSegmentAtomCount);
        Vmny.resize(maxSegmentAtomCount);
        Vmnz.resize(maxSegmentAtomCount);
        Vmn.resize(maxSegmentAtomCount * 3);
        for (int i = 0; i < maxSegmentAtomCount; i++) {
            Vmnx[i] = std::make_shared<Matrix>("Vmnx", nso, nso);
            Vmny[i] = std::make_shared<Matrix>("Vmny", nso, nso);
            Vmnz[i] = std::make_shared<Matrix>("Vmnz", nso, nso);
            Vmn[i * 3 + 0] = Vmnx[i]->get_pointer();
            Vmn[i * 3 + 1] = Vmny[i]->get_pointer();
            Vmn[i * 3 + 2] = Vmnz[i]->get_pointer();
        }
    }
#endif

    const auto& shell_pairs = Vint->shellpairs();
    size_t n_pairs = shell_pairs.size();
    for (int A = 0; A < natom; A++) {
#ifdef USING_BrianQC
        if (brianEnable) {
            if (segmentAtomCount < 0 || A < segmentAtomIndexStart || A >= (segmentAtomIndexStart + segmentAtomCount)) {
                segmentAtomIndexStart = A;
                segmentAtomCount = (segmentAtomIndexStart + maxSegmentAtomCount > natom) ? (natom - segmentAtomIndexStart) : maxSegmentAtomCount;
                
                brianInt integralType = BRIAN_INTEGRAL_TYPE_NUCLEAR;
                brianCPHFBuild1eDeriv(&brianCookie, &integralType, &segmentAtomCount, &segmentAtomIndexStart, Vmn.data());
            }
            
            C_DGEMM('N', 'N', nso, nocc, nso, 1.0, Vmnx[A - segmentAtomIndexStart]->get_pointer(), nso, Cocc->get_pointer(), nocc, 0.0, Vmix->get_pointer(), nocc);
            C_DGEMM('N', 'N', nso, nocc, nso, 1.0, Vmny[A - segmentAtomIndexStart]->get_pointer(), nso, Cocc->get_pointer(), nocc, 0.0, Vmiy->get_pointer(), nocc);
            C_DGEMM('N', 'N', nso, nocc, nso, 1.0, Vmnz[A - segmentAtomIndexStart]->get_pointer(), nso, Cocc->get_pointer(), nocc, 0.0, Vmiz->get_pointer(), nocc);
        } else {
#endif
        Vmix->zero();
        Vmiy->zero();
        Vmiz->zero();

        for (size_t p = 0; p < n_pairs; ++p) {
            auto P = shell_pairs[p].first;
            auto Q = shell_pairs[p].second;
            const auto & shellP = basisset_->shell(P);
            const auto & shellQ = basisset_->shell(Q);
            int aP = shellP.ncenter();
            int aQ = shellQ.ncenter();
            Vint->compute_shell_deriv1(P, Q);
            const auto &buffers = Vint->buffers();
            int nP = shellP.nfunction();
            int nQ = shellQ.nfunction();
            int oP = shellP.function_index();
            int oQ = shellQ.function_index();

            // buffer ordering is [Px, Py, Pz, Qx, Qy, Qz, A1x, A1y, A1z, A1x, ... ANy, ANz]
            // where the Px is the x derivative of shell P and A1x is the derivative w.r.t.
            // the nuclear charge located on atom1.  There are therefore 6 + 3*natoms buffers.
            const double* buf_x = buffers[3 * A + 6];
            const double* buf_y = buffers[3 * A + 7];
            const double* buf_z = buffers[3 * A + 8];
            double scale = P == Q ? 0.5 : 1.0;
            // Ax
            for (int p = 0; p < nP; p++) {
                for (int q = 0; q < nQ; q++) {
                    C_DAXPY(nocc, scale * (*buf_x), Cop[q + oQ], 1, Vmixp[p + oP], 1);
                    C_DAXPY(nocc, scale * (*buf_x++), Cop[p + oP], 1, Vmixp[q + oQ], 1);
                }
            }

            // Ay
            for (int p = 0; p < nP; p++) {
                for (int q = 0; q < nQ; q++) {
                    C_DAXPY(nocc, scale * (*buf_y), Cop[q + oQ], 1, Vmiyp[p + oP], 1);
                    C_DAXPY(nocc, scale * (*buf_y++), Cop[p + oP], 1, Vmiyp[q + oQ], 1);
                }
            }

            // Az
            for (int p = 0; p < nP; p++) {
                for (int q = 0; q < nQ; q++) {
                    C_DAXPY(nocc, scale * (*buf_z), Cop[q + oQ], 1, Vmizp[p + oP], 1);
                    C_DAXPY(nocc, scale * (*buf_z++), Cop[p + oP], 1, Vmizp[q + oQ], 1);
                }
            }

            // (P| derivatives
            if (aP == A) {
                const double* buf_x = buffers[0];
                const double* buf_y = buffers[1];
                const double* buf_z = buffers[2];
                // Px
                for (int p = 0; p < nP; p++) {
                    for (int q = 0; q < nQ; q++) {
                        C_DAXPY(nocc, scale * (*buf_x), Cop[q + oQ], 1, Vmixp[p + oP], 1);
                        C_DAXPY(nocc, scale * (*buf_x++), Cop[p + oP], 1, Vmixp[q + oQ], 1);
                    }
                }

                // Py
                for (int p = 0; p < nP; p++) {
                    for (int q = 0; q < nQ; q++) {
                        C_DAXPY(nocc, scale * (*buf_y), Cop[q + oQ], 1, Vmiyp[p + oP], 1);
                        C_DAXPY(nocc, scale * (*buf_y++), Cop[p + oP], 1, Vmiyp[q + oQ], 1);
                    }
                }

                // Pz
                for (int p = 0; p < nP; p++) {
                    for (int q = 0; q < nQ; q++) {
                        C_DAXPY(nocc, scale * (*buf_z), Cop[q + oQ], 1, Vmizp[p + oP], 1);
                        C_DAXPY(nocc, scale * (*buf_z++), Cop[p + oP], 1, Vmizp[q + oQ], 1);
                    }
                }
            }

            // |Q) derivatives
            if (aQ == A) {
                const double* buf_x = buffers[3];
                const double* buf_y = buffers[4];
                const double* buf_z = buffers[5];
                // Qx
                for (int p = 0; p < nP; p++) {
                    for (int q = 0; q < nQ; q++) {
                        C_DAXPY(nocc, scale * (*buf_x), Cop[q + oQ], 1, Vmixp[p + oP], 1);
                        C_DAXPY(nocc, scale * (*buf_x++), Cop[p + oP], 1, Vmixp[q + oQ], 1);
                    }
                }

                // Qy
                for (int p = 0; p < nP; p++) {
                    for (int q = 0; q < nQ; q++) {
                        C_DAXPY(nocc, scale * (*buf_y), Cop[q + oQ], 1, Vmiyp[p + oP], 1);
                        C_DAXPY(nocc, scale * (*buf_y++), Cop[p + oP], 1, Vmiyp[q + oQ], 1);
                    }
                }

                // Qz
                for (int p = 0; p < nP; p++) {
                    for (int q = 0; q < nQ; q++) {
                        C_DAXPY(nocc, scale * (*buf_z), Cop[q + oQ], 1, Vmizp[p + oP], 1);
                        C_DAXPY(nocc, scale * (*buf_z++), Cop[p + oP], 1, Vmizp[q + oQ], 1);
                    }
                }
            }
        }
#ifdef USING_BrianQC
        }
#endif
        auto Vpi_str = (alpha) ? "Vpi^A_a" : "Vpi^A_b";
        // Vpi_x
        C_DGEMM('T','N',nmo,nocc,nso,1.0,Cp[0],nmo,Vmixp[0],nocc,0.0,Vpip[0],nocc);
        psio_->write(PSIF_HESS,Vpi_str,(char*)Vpip[0], static_cast<size_t> (nmo) * nocc * sizeof(double),next_Vpi,&next_Vpi);
        // Vpi_y
        C_DGEMM('T','N',nmo,nocc,nso,1.0,Cp[0],nmo,Vmiyp[0],nocc,0.0,Vpip[0],nocc);
        psio_->write(PSIF_HESS,Vpi_str,(char*)Vpip[0], static_cast<size_t> (nmo) * nocc * sizeof(double),next_Vpi,&next_Vpi);
        // Vpi_z
        C_DGEMM('T','N',nmo,nocc,nso,1.0,Cp[0],nmo,Vmizp[0],nocc,0.0,Vpip[0],nocc);
        psio_->write(PSIF_HESS,Vpi_str,(char*)Vpip[0], static_cast<size_t> (nmo) * nocc * sizeof(double),next_Vpi,&next_Vpi);
    }
}

void USCFDeriv::JK_deriv1(std::shared_ptr<Matrix> D1,
                          std::shared_ptr<Matrix> C1, 
                          std::shared_ptr<Matrix> C1occ,
                          std::shared_ptr<Matrix> D2, 
                          int nso, int nocc, int nvir, bool alpha)
{
    // Computed JK deriv. for spin 1
    double** D1p  = D1->pointer();  
    double** C1p  = C1->pointer();  
    double** C1op = C1occ->pointer();

    double** D2p  = D2->pointer();  

    auto Dt = D1->clone();
    Dt->add(D2);
    double** Dtp = Dt->pointer();

    size_t nmo = nocc + nvir;
    int natom = molecule_->natom();

    // Figure out DFT functional info
    double Kscale = functional_->x_alpha();
    if (functional_->is_x_lrc())
        throw PSIEXCEPTION("Hessians for LRC functionals are not implemented yet.");

    size_t memory = 0.9 * memory_ / 8L;
    size_t max_a = memory / (3L * nso * nso);
    max_a = (max_a > 3 * natom ? 3 * natom : max_a);

    std::vector<SharedMatrix> dGmats;
    std::vector<double**> pdG(3*natom);
    std::vector<bool> pert_incore(3*natom);
    for(int a = 0; a < max_a; ++a)
        dGmats.push_back(std::make_shared<Matrix>("G derivative contribution", nso, nso));

    auto Gpi = std::make_shared<Matrix>("MO G Deriv", nmo, nocc);
    double** pGpi = Gpi->pointer();
    psio_address next_Gpi = PSIO_ZERO;
    auto Gpi_str = (alpha) ? "Gpi^A_a" : "Gpi^A_b";

    // Write some junk for now, to set the sizing for PSIO
    for(int pert = 0; pert < 3*natom; ++pert){
        psio_->write(PSIF_HESS,Gpi_str,(char*)pGpi[0], static_cast<size_t> (nmo) * nocc * sizeof(double),next_Gpi,&next_Gpi);
    }

    next_Gpi = PSIO_ZERO;

    if (options_.get_str("SCF_TYPE").find("DF") != std::string::npos){
        /*
         *  The DF algorithm
         */
        std::shared_ptr<BasisSet> auxiliary_ = get_basisset("DF_BASIS_SCF");

        auto Pmnfactory = std::make_shared<IntegralFactory>(auxiliary_, BasisSet::zero_ao_basis_set(), basisset_, basisset_);
        auto PQfactory = std::make_shared<IntegralFactory>(auxiliary_, BasisSet::zero_ao_basis_set(), auxiliary_, BasisSet::zero_ao_basis_set());
        std::shared_ptr<TwoBodyAOInt> Pmnint(Pmnfactory->eri(2));
        std::shared_ptr<TwoBodyAOInt> PQint(PQfactory->eri(2));
        int np = auxiliary_->nbf();
        int nso = basisset_->nbf();
        int nauxshell = auxiliary_->nshell();
        int nshell = basisset_->nshell();
        int maxp = auxiliary_->max_function_per_shell();

        auto Amn = std::make_shared<Matrix>("(A|mn)", np, nso*nso);
        auto Ami = std::make_shared<Matrix>("(A|mi)", np, nso*nocc);
        auto Aij = std::make_shared<Matrix>("(A|ij)", np, nocc*nocc);
        auto Bmn = std::make_shared<Matrix>("Minv[B][A] (A|mn)", np, nso*nso);
        auto Tmn = std::make_shared<Matrix>("Tmn", np, nso*nso);
        auto TempP = std::make_shared<Matrix>("Temp[P]", 9, maxp);
        auto TempPmn = std::make_shared<Matrix>("Temp[P][mn]", maxp, nso*nso);
        auto c = std::make_shared<Vector>("c[A] = (mn|A) D[m][n]", np);
        auto d = std::make_shared<Vector>("d[A] = Minv[A][B] C[B]", np);
        double **Amnp = Amn->pointer();
        double **Amip = Ami->pointer();
        double **Aijp = Aij->pointer();
        double **Bmnp = Bmn->pointer();
        double **pTmn = Tmn->pointer();
        double **pTempP = TempP->pointer();
        double **pTmpPmn = TempPmn->pointer();
        double *cp = c->pointer();
        double *dp = d->pointer();

        // This probably shouldn't be recomputed here; we already needed it to get the
        // second derivative integrals.  One fine day, this should be refactored.
        auto metric = std::make_shared<FittingMetric>(auxiliary_, true);            
        metric->form_full_eig_inverse(options_.get_double("DF_FITTING_CONDITION"));
        SharedMatrix PQ = metric->get_metric();
        double** PQp = PQ->pointer();

        // Same applies to these terms.  There are already hooks to compute c and d vectors, and store them on disk,
        // so we should make better use of those intermediates between the second derivative integrals and these
        // first derivative terms needed for the Fock matrix derivatives.
        for (int P = 0; P < nauxshell; ++P){
            int nP = auxiliary_->shell(P).nfunction();
            int oP = auxiliary_->shell(P).function_index();
            for(int M = 0; M < nshell; ++M){
                int nM = basisset_->shell(M).nfunction();
                int oM = basisset_->shell(M).function_index();
                for(int N = 0; N < nshell; ++N){
                    int nN = basisset_->shell(N).nfunction();
                    int oN = basisset_->shell(N).function_index();

                    Pmnint->compute_shell(P,0,M,N);
                    const double* buffer = Pmnint->buffers()[0];

                    for (int p = oP; p < oP+nP; p++) {
                        for (int m = oM; m < oM+nM; m++) {
                            for (int n = oN; n < oN+nN; n++) {
                                Amnp[p][m*nso+n] = (*buffer++);
                            }
                        }
                    }
                }
            }
        }
        // c[A] = (A|mn) D[m][n]
        C_DGEMV('N', np, nso*(size_t)nso, 1.0, Amnp[0], nso*(size_t)nso, Dtp[0], 1, 0.0, cp, 1);
        // d[A] = Minv[A][B] c[B]
        //C_DGEMV('n', np, np, 2.0, PQp[0], np, cp, 1, 0.0, dp, 1);
        C_DGEMV('n', np, np, 1.0, PQp[0], np, cp, 1, 0.0, dp, 1);

        // B[B][m,n] = Minv[A][B] (A|mn)
        C_DGEMM('n','n', np, nso*nso, np, 1.0, PQp[0], np, Amnp[0], nso*nso, 0.0, Bmnp[0], nso*nso);

        // This intermediate is for K, so use D1 rather than Dt
        // T[p][m,n] = B[p][r,n] D[m,r]
#pragma omp parallel for
        for(int p = 0; p < np; ++p)
            C_DGEMM('t', 'n', nso, nso, nso, 1.0, D1p[0], nso, Bmnp[p], nso, 0.0, pTmn[p], nso);


        for (int A = 0; A < 3 * natom; A+=max_a) {
            int nA = (A + max_a >= 3 * natom ? 3 * natom - A : max_a);

            // Keep track of which centers are loaded into memory, so we know when to skip
            std::fill(pert_incore.begin(), pert_incore.end(), false);
            std::fill(pdG.begin(), pdG.end(), (double**)nullptr);
            for (int a = 0; a < nA; a++){
                pert_incore[floor((A+a)/3.0)] = true;
                pdG[A+a] = dGmats[a]->pointer();
                dGmats[a]->zero();
            }

            for (int P = 0; P < nauxshell; ++P){
                int nP = auxiliary_->shell(P).nfunction();
                int oP = auxiliary_->shell(P).function_index();
                int Pcenter = auxiliary_->shell(P).ncenter();
                int Pncart = auxiliary_->shell(P).ncartesian();
                int Px = 3 * Pcenter + 0;
                int Py = 3 * Pcenter + 1;
                int Pz = 3 * Pcenter + 2;
                for(int Q = 0; Q < nauxshell; ++Q){
                    int nQ = auxiliary_->shell(Q).nfunction();
                    int oQ = auxiliary_->shell(Q).function_index();
                    int Qcenter = auxiliary_->shell(Q).ncenter();
                    int Qncart = auxiliary_->shell(Q).ncartesian();
                    int Qx = 3 * Qcenter + 0;
                    int Qy = 3 * Qcenter + 1;
                    int Qz = 3 * Qcenter + 2;

                    if(!pert_incore[Pcenter] && !pert_incore[Qcenter])
                        continue;

                    PQint->compute_shell_deriv1(P,0,Q,0);
                    const auto& buffers = PQint->buffers();
                    const double* PxBuf = buffers[0];
                    const double* PyBuf = buffers[1];
                    const double* PzBuf = buffers[2];
                    const double* QxBuf = buffers[3];
                    const double* QyBuf = buffers[4];
                    const double* QzBuf = buffers[5];

                    if(pert_incore[Pcenter]){
                        // J terms
                        // Px
                        C_DGEMV('n', nP, nQ, 1.0, const_cast<double*>(PxBuf), nQ, &dp[oQ], 1, 0.0, pTempP[0], 1);
                        C_DGEMV('t', nP, nso*nso, -1.0, Bmnp[oP], nso*nso, pTempP[0], 1, 1.0, pdG[Px][0], 1);
                        // Py
                        C_DGEMV('n', nP, nQ, 1.0, const_cast<double*>(PyBuf), nQ, &dp[oQ], 1, 0.0, pTempP[0], 1);
                        C_DGEMV('t', nP, nso*nso, -1.0, Bmnp[oP], nso*nso, pTempP[0], 1, 1.0, pdG[Py][0], 1);
                        // Pz
                        C_DGEMV('n', nP, nQ, 1.0, const_cast<double*>(PzBuf), nQ, &dp[oQ], 1, 0.0, pTempP[0], 1);
                        C_DGEMV('t', nP, nso*nso, -1.0, Bmnp[oP], nso*nso, pTempP[0], 1, 1.0, pdG[Pz][0], 1);

                        if (Kscale) {
                            // K terms
                            // Px
                            C_DGEMM('n', 'n', nP, nso*nso, nQ, 1.0, const_cast<double*>(PxBuf), nQ, pTmn[oQ], nso*nso, 0.0, pTmpPmn[0], nso*nso);
                            for(int p = 0; p < nP; ++p)
                                C_DGEMM('N', 'N', nso, nso, nso, Kscale, Bmnp[p+oP], nso, pTmpPmn[p], nso, 1.0, pdG[Px][0], nso);
                            // Py
                            C_DGEMM('n', 'n', nP, nso*nso, nQ, 1.0, const_cast<double*>(PyBuf), nQ, pTmn[oQ], nso*nso, 0.0, pTmpPmn[0], nso*nso);
                            for(int p = 0; p < nP; ++p)
                                C_DGEMM('N', 'N', nso, nso, nso, Kscale, Bmnp[p+oP], nso, pTmpPmn[p], nso, 1.0, pdG[Py][0], nso);
                            // Pz
                            C_DGEMM('n', 'n', nP, nso*nso, nQ, 1.0, const_cast<double*>(PzBuf), nQ, pTmn[oQ], nso*nso, 0.0, pTmpPmn[0], nso*nso);
                            for(int p = 0; p < nP; ++p)
                                C_DGEMM('N', 'N', nso, nso, nso, Kscale, Bmnp[p+oP], nso, pTmpPmn[p], nso, 1.0, pdG[Pz][0], nso);
                        }

                    }
                    if(pert_incore[Qcenter]){
                        // J terms
                        // Qx
                        C_DGEMV('n', nP, nQ, 1.0, const_cast<double*>(QxBuf), nQ, &dp[oQ], 1, 0.0, pTempP[0], 1);
                        C_DGEMV('t', nP, nso*nso, -1.0, Bmnp[oP], nso*nso, pTempP[0], 1, 1.0, pdG[Qx][0], 1);
                        // Qy
                        C_DGEMV('n', nP, nQ, 1.0, const_cast<double*>(QyBuf), nQ, &dp[oQ], 1, 0.0, pTempP[0], 1);
                        C_DGEMV('t', nP, nso*nso, -1.0, Bmnp[oP], nso*nso, pTempP[0], 1, 1.0, pdG[Qy][0], 1);
                        // Qz
                        C_DGEMV('n', nP, nQ, 1.0, const_cast<double*>(QzBuf), nQ, &dp[oQ], 1, 0.0, pTempP[0], 1);
                        C_DGEMV('t', nP, nso*nso, -1.0, Bmnp[oP], nso*nso, pTempP[0], 1, 1.0, pdG[Qz][0], 1);

                        if (Kscale) {
                            // K terms
                            // Qx
                            C_DGEMM('n', 'n', nP, nso*nso, nQ, 1.0, const_cast<double*>(QxBuf), nQ, pTmn[oQ], nso*nso, 0.0, pTmpPmn[0], nso*nso);
                            for(int p = 0; p < nP; ++p)
                                C_DGEMM('N', 'N', nso, nso, nso, Kscale, Bmnp[p+oP], nso, pTmpPmn[p], nso, 1.0, pdG[Qx][0], nso);
                            // Qy
                            C_DGEMM('n', 'n', nP, nso*nso, nQ, 1.0, const_cast<double*>(QyBuf), nQ, pTmn[oQ], nso*nso, 0.0, pTmpPmn[0], nso*nso);
                            for(int p = 0; p < nP; ++p)
                                C_DGEMM('N', 'N', nso, nso, nso, Kscale, Bmnp[p+oP], nso, pTmpPmn[p], nso, 1.0, pdG[Qy][0], nso);
                            // Qz
                            C_DGEMM('n', 'n', nP, nso*nso, nQ, 1.0, const_cast<double*>(QzBuf), nQ, pTmn[oQ], nso*nso, 0.0, pTmpPmn[0], nso*nso);
                            for(int p = 0; p < nP; ++p)
                                C_DGEMM('N', 'N', nso, nso, nso, Kscale, Bmnp[p+oP], nso, pTmpPmn[p], nso, 1.0, pdG[Qz][0], nso);
                        }
                    }

                }
            }


            for (int P = 0; P < nauxshell; ++P){
                int nP = auxiliary_->shell(P).nfunction();
                int oP = auxiliary_->shell(P).function_index();
                int Pcenter = auxiliary_->shell(P).ncenter();
                int Pncart = auxiliary_->shell(P).ncartesian();
                int Px = 3 * Pcenter + 0;
                int Py = 3 * Pcenter + 1;
                int Pz = 3 * Pcenter + 2;
                for(int M = 0; M < nshell; ++M){
                    int nM = basisset_->shell(M).nfunction();
                    int oM = basisset_->shell(M).function_index();
                    int Mcenter = basisset_->shell(M).ncenter();
                    int Mncart = basisset_->shell(M).ncartesian();
                    int mx = 3 * Mcenter + 0;
                    int my = 3 * Mcenter + 1;
                    int mz = 3 * Mcenter + 2;
                    for(int N = 0; N < nshell; ++N){
                        int nN = basisset_->shell(N).nfunction();
                        int oN = basisset_->shell(N).function_index();
                        int Ncenter = basisset_->shell(N).ncenter();
                        int Nncart = basisset_->shell(N).ncartesian();
                        int nx = 3 * Ncenter + 0;
                        int ny = 3 * Ncenter + 1;
                        int nz = 3 * Ncenter + 2;

                        if(!pert_incore[Pcenter] &&
                                !pert_incore[Mcenter] &&
                                !pert_incore[Ncenter])
                            continue;

                        Pmnint->compute_shell_deriv1(P,0,M,N);
                        const auto& buffers = Pmnint->buffers();
                        const double* PxBuf = buffers[0];
                        const double* PyBuf = buffers[1];
                        const double* PzBuf = buffers[2];
                        const double* mxBuf = buffers[3];
                        const double* myBuf = buffers[4];
                        const double* mzBuf = buffers[5];
                        const double* nxBuf = buffers[6];
                        const double* nyBuf = buffers[7];
                        const double* nzBuf = buffers[8];

                        /*
                         * J terms have 2 contributions:
                         *      F^x[m][n] <- (P|mn)^x d[P]
                         * and
                         *      F^x[r][s] <- D[m][n] (P|mn)^x B[P][r,s]
                         * The second term factorizes into...
                         * ... Temp[P] = D[m][n] (P|mn)^x ...
                         * ... and then F^x[r][s] <- Temp[P] B[P][r,s]
                         */

                        for(int x = 0; x < 9; ++ x){
                            size_t delta = 0L;
                            for(int p = 0; p < nP; ++p){
                                double val = 0.0;
                                for(int m = oM; m < nM+oM; ++m){
                                    for(int n = oN; n < nN+oN; ++n){
                                        val += Dtp[m][n] * buffers[x][delta];
                                        ++delta;
                                    }
                                }
                                pTempP[x][p] = val;
                            }
                        }

                        if(pert_incore[Pcenter]){
                            // J Terms
                            size_t delta = 0L;
                            for(int p = oP; p < oP+nP; ++p){
                                for(int m = oM; m < nM+oM; ++m){
                                    for(int n = oN; n < nN+oN; ++n){
                                        pdG[Px][m][n] += PxBuf[delta] * dp[p];
                                        pdG[Py][m][n] += PyBuf[delta] * dp[p];
                                        pdG[Pz][m][n] += PzBuf[delta] * dp[p];
                                        ++delta;
                                    }
                                }
                            }
                            C_DGEMV('t', nP, nso*nso, 1.0, Bmnp[oP], nso*nso, pTempP[0], 1, 1.0, pdG[Px][0], 1);
                            C_DGEMV('t', nP, nso*nso, 1.0, Bmnp[oP], nso*nso, pTempP[1], 1, 1.0, pdG[Py][0], 1);
                            C_DGEMV('t', nP, nso*nso, 1.0, Bmnp[oP], nso*nso, pTempP[2], 1, 1.0, pdG[Pz][0], 1);
                            // K Terms
                            if (Kscale) {
                                for(int p = 0; p < nP; ++p)
                                    C_DGEMM('T', 'N', nN, nso, nM, -2 * Kscale, const_cast<double*>(&PxBuf[p*nM*nN]), nN, &pTmn[oP+p][oM*nso], nso, 1.0, pdG[Px][oN], nso);
                                for(int p = 0; p < nP; ++p)
                                    C_DGEMM('T', 'N', nN, nso, nM, -2 * Kscale, const_cast<double*>(&PyBuf[p*nM*nN]), nN, &pTmn[oP+p][oM*nso], nso, 1.0, pdG[Py][oN], nso);
                                for(int p = 0; p < nP; ++p)
                                    C_DGEMM('T', 'N', nN, nso, nM, -2 * Kscale, const_cast<double*>(&PzBuf[p*nM*nN]), nN, &pTmn[oP+p][oM*nso], nso, 1.0, pdG[Pz][oN], nso);
                            }
                        }
                        if(pert_incore[Mcenter]){
                            // J Terms
                            size_t delta = 0L;
                            for(int p = oP; p < oP+nP; ++p){
                                for(int m = oM; m < nM+oM; ++m){
                                    for(int n = oN; n < nN+oN; ++n){
                                        pdG[mx][m][n] += mxBuf[delta] * dp[p];
                                        pdG[my][m][n] += myBuf[delta] * dp[p];
                                        pdG[mz][m][n] += mzBuf[delta] * dp[p];
                                        ++delta;
                                    }
                                }
                            }
                            C_DGEMV('t', nP, nso*nso, 1.0, Bmnp[oP], nso*nso, pTempP[3], 1, 1.0, pdG[mx][0], 1);
                            C_DGEMV('t', nP, nso*nso, 1.0, Bmnp[oP], nso*nso, pTempP[4], 1, 1.0, pdG[my][0], 1);
                            C_DGEMV('t', nP, nso*nso, 1.0, Bmnp[oP], nso*nso, pTempP[5], 1, 1.0, pdG[mz][0], 1);
                            // K Terms
                            if (Kscale) {
                                for(int p = 0; p < nP; ++p)
                                    C_DGEMM('T', 'N', nN, nso, nM, -2 * Kscale, const_cast<double*>(&mxBuf[p*nM*nN]), nN, &pTmn[oP+p][oM*nso], nso, 1.0, pdG[mx][oN], nso);
                                for(int p = 0; p < nP; ++p)
                                    C_DGEMM('T', 'N', nN, nso, nM, -2 * Kscale, const_cast<double*>(&myBuf[p*nM*nN]), nN, &pTmn[oP+p][oM*nso], nso, 1.0, pdG[my][oN], nso);
                                for(int p = 0; p < nP; ++p)
                                    C_DGEMM('T', 'N', nN, nso, nM, -2 * Kscale, const_cast<double*>(&mzBuf[p*nM*nN]), nN, &pTmn[oP+p][oM*nso], nso, 1.0, pdG[mz][oN], nso);
                            }
                        }
                        if(pert_incore[Ncenter]){
                            // J Terms
                            size_t delta = 0L;
                            for(int p = oP; p < oP+nP; ++p){
                                for(int m = oM; m < nM+oM; ++m){
                                    for(int n = oN; n < nN+oN; ++n){
                                        pdG[nx][m][n] += nxBuf[delta] * dp[p];
                                        pdG[ny][m][n] += nyBuf[delta] * dp[p];
                                        pdG[nz][m][n] += nzBuf[delta] * dp[p];
                                        ++delta;
                                    }
                                }
                            }
                            C_DGEMV('t', nP, nso*nso, 1.0, Bmnp[oP], nso*nso, pTempP[6], 1, 1.0, pdG[nx][0], 1);
                            C_DGEMV('t', nP, nso*nso, 1.0, Bmnp[oP], nso*nso, pTempP[7], 1, 1.0, pdG[ny][0], 1);
                            C_DGEMV('t', nP, nso*nso, 1.0, Bmnp[oP], nso*nso, pTempP[8], 1, 1.0, pdG[nz][0], 1);
                            // K Terms
                            if (Kscale) {
                                for(int p = 0; p < nP; ++p)
                                    C_DGEMM('T', 'N', nN, nso, nM, -2 * Kscale, const_cast<double*>(&nxBuf[p*nM*nN]), nN, &pTmn[oP+p][oM*nso], nso, 1.0, pdG[nx][oN], nso);
                                for(int p = 0; p < nP; ++p)
                                    C_DGEMM('T', 'N', nN, nso, nM, -2 * Kscale, const_cast<double*>(&nyBuf[p*nM*nN]), nN, &pTmn[oP+p][oM*nso], nso, 1.0, pdG[ny][oN], nso);
                                for(int p = 0; p < nP; ++p)
                                    C_DGEMM('T', 'N', nN, nso, nM, -2 * Kscale, const_cast<double*>(&nzBuf[p*nM*nN]), nN, &pTmn[oP+p][oM*nso], nso, 1.0, pdG[nz][oN], nso);
                            }
                        }

                    }
                }
            }

            for(int a = 0; a < nA; ++a){
                // Symmetrize the derivative Fock contributions
                SharedMatrix G = dGmats[a];
                G->add(G->transpose());
                Gpi->transform(C1, G, C1occ);
                Gpi->scale(0.5);
                psio_->write(PSIF_HESS,Gpi_str,(char*)pGpi[0], static_cast<size_t> (nmo) * nocc * sizeof(double),next_Gpi,&next_Gpi);
            }
        } // End loop over A batches

    }else{
        /*
         * The conventional integral algorithm
         */
#ifdef USING_BrianQC
        if (brianEnable) {
            brianBool computeCoulomb = BRIAN_TRUE;
            brianBool computeExchange = BRIAN_TRUE;
            
            brianInt maxSegmentSize;
            brianCPHFMaxSegmentSize(&brianCookie, &maxSegmentSize);
            brianInt maxSegmentAtomCount = maxSegmentSize / 3;
            brianInt segmentAtomCount = -1;
            brianInt segmentAtomIndexStart = -1;
            
            std::vector<std::shared_ptr<Matrix>> Jmnx(maxSegmentAtomCount);
            std::vector<std::shared_ptr<Matrix>> Jmny(maxSegmentAtomCount);
            std::vector<std::shared_ptr<Matrix>> Jmnz(maxSegmentAtomCount);
            std::vector<std::shared_ptr<Matrix>> Kmnx1(maxSegmentAtomCount);
            std::vector<std::shared_ptr<Matrix>> Kmny1(maxSegmentAtomCount);
            std::vector<std::shared_ptr<Matrix>> Kmnz1(maxSegmentAtomCount);
            std::vector<std::shared_ptr<Matrix>> Kmnx2(maxSegmentAtomCount);
            std::vector<std::shared_ptr<Matrix>> Kmny2(maxSegmentAtomCount);
            std::vector<std::shared_ptr<Matrix>> Kmnz2(maxSegmentAtomCount);
            std::vector<double*> Jmn(maxSegmentAtomCount * 3);
            std::vector<double*> Kmn1(maxSegmentAtomCount * 3);
            std::vector<double*> Kmn2(maxSegmentAtomCount * 3);
            for (int i = 0; i < maxSegmentAtomCount; i++) {
                Jmnx[i] = std::make_shared<Matrix>("Jmnx", nso, nso);
                Jmny[i] = std::make_shared<Matrix>("Jmny", nso, nso);
                Jmnz[i] = std::make_shared<Matrix>("Jmnz", nso, nso);
                Kmnx1[i] = std::make_shared<Matrix>("Kmnx1", nso, nso);
                Kmny1[i] = std::make_shared<Matrix>("Kmny1", nso, nso);
                Kmnz1[i] = std::make_shared<Matrix>("Kmnz1", nso, nso);
                Kmnx2[i] = std::make_shared<Matrix>("Kmnx2", nso, nso);
                Kmny2[i] = std::make_shared<Matrix>("Kmny2", nso, nso);
                Kmnz2[i] = std::make_shared<Matrix>("Kmnz2", nso, nso);
                Jmn[i * 3 + 0] = Jmnx[i]->get_pointer();
                Jmn[i * 3 + 1] = Jmny[i]->get_pointer();
                Jmn[i * 3 + 2] = Jmnz[i]->get_pointer();
                Kmn1[i * 3 + 0] = Kmnx1[i]->get_pointer();
                Kmn1[i * 3 + 1] = Kmny1[i]->get_pointer();
                Kmn1[i * 3 + 2] = Kmnz1[i]->get_pointer();
                Kmn2[i * 3 + 0] = Kmnx2[i]->get_pointer();
                Kmn2[i * 3 + 1] = Kmny2[i]->get_pointer();
                Kmn2[i * 3 + 2] = Kmnz2[i]->get_pointer();
            }
            
            for (int A = 0; A < natom; A++) {
                if (segmentAtomCount < 0 || A < segmentAtomIndexStart || A >= (segmentAtomIndexStart + segmentAtomCount)) {
                    segmentAtomIndexStart = A;
                    segmentAtomCount = (segmentAtomIndexStart + maxSegmentAtomCount > natom) ? (natom - segmentAtomIndexStart) : maxSegmentAtomCount;
                    
                    brianInt integralType = BRIAN_INTEGRAL_TYPE_NUCLEAR;
                    brianCPHFBuildRepulsionDeriv(&brianCookie, &computeCoulomb, &computeExchange, &segmentAtomCount, &segmentAtomIndexStart, D1->get_pointer(), D2->get_pointer(), Jmn.data(), Kmn1.data(), Kmn2.data());
                }
                
                Jmnx[A - segmentAtomIndexStart]->subtract(Kmnx1[A - segmentAtomIndexStart]);
                Jmny[A - segmentAtomIndexStart]->subtract(Kmny1[A - segmentAtomIndexStart]);
                Jmnz[A - segmentAtomIndexStart]->subtract(Kmnz1[A - segmentAtomIndexStart]);
                
                Gpi->transform(C1, Jmnx[A - segmentAtomIndexStart], C1occ);
                psio_->write(PSIF_HESS,Gpi_str,(char*)pGpi[0],nmo * nocc * sizeof(double),next_Gpi,&next_Gpi);
                Gpi->transform(C1, Jmny[A - segmentAtomIndexStart], C1occ);
                psio_->write(PSIF_HESS,Gpi_str,(char*)pGpi[0],nmo * nocc * sizeof(double),next_Gpi,&next_Gpi);
                Gpi->transform(C1, Jmnz[A - segmentAtomIndexStart], C1occ);
                psio_->write(PSIF_HESS,Gpi_str,(char*)pGpi[0],nmo * nocc * sizeof(double),next_Gpi,&next_Gpi);
            }
        } else {
#endif
        std::shared_ptr<TwoBodyAOInt> ints(integral_->eri(1));

        const std::vector<std::pair<int, int> >& shell_pairs = ints->shell_pairs();
        size_t npairs = shell_pairs.size();
        size_t npairs2 = npairs * npairs;

        for (int A = 0; A < 3 * natom; A+=max_a) {
            int nA = (A + max_a >= 3 * natom ? 3 * natom - A : max_a);

            // Keep track of which centers are loaded into memory, so we know when to skip
            std::fill(pert_incore.begin(), pert_incore.end(), false);
            std::fill(pdG.begin(), pdG.end(), (double**)nullptr);
            for (int a = 0; a < nA; a++){
                pert_incore[floor((A+a)/3.0)] = true;
                pdG[A+a] = dGmats[a]->pointer();
                dGmats[a]->zero();
            }

            for (size_t index = 0L; index < npairs2; index++) {

                size_t PQ = index / npairs;
                size_t RS = index % npairs;

                if (RS > PQ) continue;

                int P = shell_pairs[PQ].first;
                int Q = shell_pairs[PQ].second;
                int R = shell_pairs[RS].first;
                int S = shell_pairs[RS].second;

                if (!ints->shell_significant(P,Q,R,S)) continue;

                int Pcenter = basisset_->shell(P).ncenter();
                int Qcenter = basisset_->shell(Q).ncenter();
                int Rcenter = basisset_->shell(R).ncenter();
                int Scenter = basisset_->shell(S).ncenter();
                if(!pert_incore[Pcenter] &&
                        !pert_incore[Qcenter] &&
                        !pert_incore[Rcenter] &&
                        !pert_incore[Scenter])
                    continue;

                int am1 = basisset_->shell(P).am();
                int am2 = basisset_->shell(Q).am();
                int am3 = basisset_->shell(R).am();
                int am4 = basisset_->shell(S).am();

                // Correct libint ordering to libint's expected P>=Q, R>=S, PQ<=RS, for efficiency.
                if (am1 < am2){
                    std::swap(P,Q);
                    std::swap(Pcenter,Qcenter);
                }
                if (am3 < am4){
                    std::swap(R,S);
                    std::swap(Rcenter,Scenter);
                }
                if( (am1 + am2) > (am3 + am4) ){
                    std::swap(P,R);
                    std::swap(Q,S);
                    std::swap(Pcenter,Rcenter);
                    std::swap(Qcenter,Scenter);
                }

                ints->compute_shell_deriv1(P,Q,R,S);
                const auto& buffers = ints->buffers();
                const double* AxBuf = buffers[0];
                const double* AyBuf = buffers[1];
                const double* AzBuf = buffers[2];
                const double* BxBuf = buffers[3];
                const double* ByBuf = buffers[4];
                const double* BzBuf = buffers[5];
                const double* CxBuf = buffers[6];
                const double* CyBuf = buffers[7];
                const double* CzBuf = buffers[8];
                const double* DxBuf = buffers[9];
                const double* DyBuf = buffers[10];
                const double* DzBuf = buffers[11];

                const int Psize = basisset_->shell(P).nfunction();
                const int Qsize = basisset_->shell(Q).nfunction();
                const int Rsize = basisset_->shell(R).nfunction();
                const int Ssize = basisset_->shell(S).nfunction();

                const int Pncart = basisset_->shell(P).ncartesian();
                const int Qncart = basisset_->shell(Q).ncartesian();
                const int Rncart = basisset_->shell(R).ncartesian();
                const int Sncart = basisset_->shell(S).ncartesian();

                const int Poff = basisset_->shell(P).function_index();
                const int Qoff = basisset_->shell(Q).function_index();
                const int Roff = basisset_->shell(R).function_index();
                const int Soff = basisset_->shell(S).function_index();
        
                double prefactor = 1.0;
                if (P != Q)   prefactor *= 2.0;
                if (R != S)   prefactor *= 2.0;
                if (PQ != RS) prefactor *= 2.0;

                double Dpq, Drs, Dpr, Dqs, Dps, Dqr;
                size_t delta;
                double Ax, Ay, Az;
                double Bx, By, Bz;
                double Cx, Cy, Cz;
                double Dx, Dy, Dz;

                // => Coulomb Term <= //
                prefactor *= 0.5;

                // (p1q1|r1s1)
                Ax = 0.0; Ay = 0.0; Az = 0.0;
                Bx = 0.0; By = 0.0; Bz = 0.0;
                Cx = 0.0; Cy = 0.0; Cz = 0.0;
                Dx = 0.0; Dy = 0.0; Dz = 0.0;
                delta = 0L;
                for (int p = Poff; p < Poff+Psize; p++) {
                    for (int q = Qoff; q < Qoff+Qsize; q++) {
                        for (int r = Roff; r < Roff+Rsize; r++) {
                            for (int s = Soff; s < Soff+Ssize; s++) {
                                Dpq = Dtp[p][q];
                                Drs = Dtp[r][s];
                                Ax = prefactor * AxBuf[delta];
                                Ay = prefactor * AyBuf[delta];
                                Az = prefactor * AzBuf[delta];
                                Bx = prefactor * BxBuf[delta];
                                By = prefactor * ByBuf[delta];
                                Bz = prefactor * BzBuf[delta];
                                Cx = prefactor * CxBuf[delta];
                                Cy = prefactor * CyBuf[delta];
                                Cz = prefactor * CzBuf[delta];
                                Dx = prefactor * DxBuf[delta];
                                Dy = prefactor * DyBuf[delta];
                                Dz = prefactor * DzBuf[delta];

                                if(pert_incore[Pcenter]){
                                    pdG[Pcenter*3+0][p][q] += Ax * Drs;
                                    pdG[Pcenter*3+0][r][s] += Ax * Dpq;
                                    pdG[Pcenter*3+1][p][q] += Ay * Drs;
                                    pdG[Pcenter*3+1][r][s] += Ay * Dpq;
                                    pdG[Pcenter*3+2][p][q] += Az * Drs;
                                    pdG[Pcenter*3+2][r][s] += Az * Dpq;
                                }
                                if(pert_incore[Qcenter]){
                                    pdG[Qcenter*3+0][p][q] += Bx * Drs;
                                    pdG[Qcenter*3+0][r][s] += Bx * Dpq;
                                    pdG[Qcenter*3+1][p][q] += By * Drs;
                                    pdG[Qcenter*3+1][r][s] += By * Dpq;
                                    pdG[Qcenter*3+2][p][q] += Bz * Drs;
                                    pdG[Qcenter*3+2][r][s] += Bz * Dpq;
                                }
                                if(pert_incore[Rcenter]){
                                    pdG[Rcenter*3+0][p][q] += Cx * Drs;
                                    pdG[Rcenter*3+0][r][s] += Cx * Dpq;
                                    pdG[Rcenter*3+1][p][q] += Cy * Drs;
                                    pdG[Rcenter*3+1][r][s] += Cy * Dpq;
                                    pdG[Rcenter*3+2][p][q] += Cz * Drs;
                                    pdG[Rcenter*3+2][r][s] += Cz * Dpq;
                                }
                                if(pert_incore[Scenter]){
                                    pdG[Scenter*3+0][p][q] += Dx * Drs;
                                    pdG[Scenter*3+0][r][s] += Dx * Dpq;
                                    pdG[Scenter*3+1][p][q] += Dy * Drs;
                                    pdG[Scenter*3+1][r][s] += Dy * Dpq;
                                    pdG[Scenter*3+2][p][q] += Dz * Drs;
                                    pdG[Scenter*3+2][r][s] += Dz * Dpq;
                                }
                                delta++;
                            }
                        }
                    }
                }

                prefactor *= 2.0;
                // => Exchange Term <= //
                if(Kscale) {
                    Ax = 0.0; Ay = 0.0; Az = 0.0;
                    Bx = 0.0; By = 0.0; Bz = 0.0;
                    Cx = 0.0; Cy = 0.0; Cz = 0.0;
                    Dx = 0.0; Dy = 0.0; Dz = 0.0;
                    delta = 0L;
                    prefactor *= -0.25 * Kscale;
                    for (int p = Poff; p < Poff+Psize; p++) {
                        for (int q = Qoff; q < Qoff+Qsize; q++) {
                            for (int r = Roff; r < Roff+Rsize; r++) {
                                for (int s = Soff; s < Soff+Ssize; s++) {
                                    Ax = prefactor * AxBuf[delta];
                                    Ay = prefactor * AyBuf[delta];
                                    Az = prefactor * AzBuf[delta];
                                    Bx = prefactor * BxBuf[delta];
                                    By = prefactor * ByBuf[delta];
                                    Bz = prefactor * BzBuf[delta];
                                    Cx = prefactor * CxBuf[delta];
                                    Cy = prefactor * CyBuf[delta];
                                    Cz = prefactor * CzBuf[delta];
                                    Dx = prefactor * DxBuf[delta];
                                    Dy = prefactor * DyBuf[delta];
                                    Dz = prefactor * DzBuf[delta];

                                    Dpr = D1p[p][r];
                                    Dqs = D1p[q][s];
                                    Dps = D1p[p][s];
                                    Dqr = D1p[q][r];
                                    if(pert_incore[Pcenter]){
                                        pdG[Pcenter*3+0][p][r] += Ax * Dqs;
                                        pdG[Pcenter*3+0][q][s] += Ax * Dpr;
                                        pdG[Pcenter*3+0][p][s] += Ax * Dqr;
                                        pdG[Pcenter*3+0][q][r] += Ax * Dps;
                                        pdG[Pcenter*3+1][p][r] += Ay * Dqs;
                                        pdG[Pcenter*3+1][q][s] += Ay * Dpr;
                                        pdG[Pcenter*3+1][p][s] += Ay * Dqr;
                                        pdG[Pcenter*3+1][q][r] += Ay * Dps;
                                        pdG[Pcenter*3+2][p][r] += Az * Dqs;
                                        pdG[Pcenter*3+2][q][s] += Az * Dpr;
                                        pdG[Pcenter*3+2][p][s] += Az * Dqr;
                                        pdG[Pcenter*3+2][q][r] += Az * Dps;
                                    }
                                    if(pert_incore[Qcenter]){
                                        pdG[Qcenter*3+0][p][r] += Bx * Dqs;
                                        pdG[Qcenter*3+0][q][s] += Bx * Dpr;
                                        pdG[Qcenter*3+0][p][s] += Bx * Dqr;
                                        pdG[Qcenter*3+0][q][r] += Bx * Dps;
                                        pdG[Qcenter*3+1][p][r] += By * Dqs;
                                        pdG[Qcenter*3+1][q][s] += By * Dpr;
                                        pdG[Qcenter*3+1][p][s] += By * Dqr;
                                        pdG[Qcenter*3+1][q][r] += By * Dps;
                                        pdG[Qcenter*3+2][p][r] += Bz * Dqs;
                                        pdG[Qcenter*3+2][q][s] += Bz * Dpr;
                                        pdG[Qcenter*3+2][p][s] += Bz * Dqr;
                                        pdG[Qcenter*3+2][q][r] += Bz * Dps;
                                    }
                                    if(pert_incore[Rcenter]){
                                        pdG[Rcenter*3+0][p][r] += Cx * Dqs;
                                        pdG[Rcenter*3+0][q][s] += Cx * Dpr;
                                        pdG[Rcenter*3+0][p][s] += Cx * Dqr;
                                        pdG[Rcenter*3+0][q][r] += Cx * Dps;
                                        pdG[Rcenter*3+1][p][r] += Cy * Dqs;
                                        pdG[Rcenter*3+1][q][s] += Cy * Dpr;
                                        pdG[Rcenter*3+1][p][s] += Cy * Dqr;
                                        pdG[Rcenter*3+1][q][r] += Cy * Dps;
                                        pdG[Rcenter*3+2][p][r] += Cz * Dqs;
                                        pdG[Rcenter*3+2][q][s] += Cz * Dpr;
                                        pdG[Rcenter*3+2][p][s] += Cz * Dqr;
                                        pdG[Rcenter*3+2][q][r] += Cz * Dps;
                                    }
                                    if(pert_incore[Scenter]){
                                        pdG[Scenter*3+0][p][r] += Dx * Dqs;
                                        pdG[Scenter*3+0][q][s] += Dx * Dpr;
                                        pdG[Scenter*3+0][p][s] += Dx * Dqr;
                                        pdG[Scenter*3+0][q][r] += Dx * Dps;
                                        pdG[Scenter*3+1][p][r] += Dy * Dqs;
                                        pdG[Scenter*3+1][q][s] += Dy * Dpr;
                                        pdG[Scenter*3+1][p][s] += Dy * Dqr;
                                        pdG[Scenter*3+1][q][r] += Dy * Dps;
                                        pdG[Scenter*3+2][p][r] += Dz * Dqs;
                                        pdG[Scenter*3+2][q][s] += Dz * Dpr;
                                        pdG[Scenter*3+2][p][s] += Dz * Dqr;
                                        pdG[Scenter*3+2][q][r] += Dz * Dps;
                                    }
                                    delta++;
                                }
                            }
                        }
                    }
                }
            } // End shell loops

            for(int a = 0; a < nA; ++a){
                // Symmetrize the derivative Fock contributions
                SharedMatrix G = dGmats[a];
                G->add(G->transpose());
                Gpi->transform(C1, G, C1occ);
                Gpi->scale(0.5);
                psio_->write(PSIF_HESS,Gpi_str,(char*)pGpi[0], static_cast<size_t> (nmo) * nocc * sizeof(double),next_Gpi,&next_Gpi);
            }
        } // End loop over A batches
#ifdef USING_BrianQC
        }
#endif
    } // End if density fitted
}

void USCFDeriv::JK_deriv2(std::shared_ptr<JK> jk, int mem,
                          std::shared_ptr<Matrix> C1, 
                          std::shared_ptr<Matrix> C1occ,
                          std::shared_ptr<Matrix> C2, 
                          std::shared_ptr<Matrix> C2occ,
                          int nso, int n1occ, int n2occ, int n1vir, bool alpha)
{
    // => J2pi/K2pi <= //
    double** C1p  = C1->pointer();  
    double** C1op = C1occ->pointer();

    double** C2p  = C2->pointer();  
    double** C2op = C2occ->pointer();
    size_t nmo = n1occ + n1vir;
    int natom = molecule_->natom();

    size_t per_A = 3L * nso * nso + 1L * n1occ * nso;
    size_t max_A = (mem / 2L) / per_A;
    max_A = (max_A > 3 * natom ? 3 * natom : max_A);

    // Figure out DFT functional info
    double Kscale = functional_->x_alpha();
    if (functional_->is_x_lrc())
        throw PSIEXCEPTION("Hessians for LRC functionals are not implemented yet.");

    auto Sij1 = std::make_shared<Matrix>("Sij1",n1occ,n1occ);
    double** Sij1p = Sij1->pointer();
    auto Sij2 = std::make_shared<Matrix>("Sij2",n2occ,n2occ);
    double** Sij2p = Sij2->pointer();
    auto T = std::make_shared<Matrix>("T",nso,n1occ);
    double** Tp = T->pointer();
    auto U = std::make_shared<Matrix>("Tempai",nmo,n1occ);
    double** Up = U->pointer();
    auto G2pi_str = (alpha) ? "G2pi^A_a" : "G2pi^A_b";
    auto Sij_1 = (alpha) ? "Sij^A_a" : "Sij^A_b";
    auto Sij_2 = (alpha) ? "Sij^A_b" : "Sij^A_a";

    std::vector<std::shared_ptr<Matrix> >& L = jk->C_left();
    std::vector<std::shared_ptr<Matrix> >& R = jk->C_right();
    const std::vector<std::shared_ptr<Matrix> >& J = jk->J();
    const std::vector<std::shared_ptr<Matrix> >& K = jk->K();
    L.clear();
    R.clear();

    // Write some placeholder data to PSIO, to get the sizing right
    psio_address next_Gpi = PSIO_ZERO;
    for (int A = 0; A < 3 * natom; A++)
        psio_->write(PSIF_HESS,G2pi_str,(char*)Up[0], static_cast<size_t> (nmo)*n1occ*sizeof(double),next_Gpi,&next_Gpi);

    std::vector<SharedMatrix> Dx, Vx;
    for (int A = 0; A < max_A; A++) {
        // Just pass C1 quantities in; this object doesn't respect symmetry anyway
        L.push_back(C1occ);
        R.push_back(std::make_shared<Matrix>("R",nso,n1occ));
        Dx.push_back(std::make_shared<Matrix>("Dx", nso,nso));
        Vx.push_back(std::make_shared<Matrix>("Vx", nso,nso));

        L.push_back(C2occ);
        R.push_back(std::make_shared<Matrix>("R",nso,n2occ));
    }

    jk->print_header();

    // First the all spin-1 part
    for (int A = 0; A < 3 * natom; A+=max_A) {
        int nA = max_A;
        if (A + max_A >= 3 * natom) {
            nA = 3 * natom - A;
            L.resize(2*nA);
            R.resize(2*nA);
        }
        for (int a = 0; a < nA; a++) {
            psio_address next_Sij= psio_get_address(PSIO_ZERO,(A + a) * (size_t) n1occ * n1occ * sizeof(double));
            psio_->read(PSIF_HESS,Sij_1,(char*)Sij1p[0], static_cast<size_t> (n1occ)*n1occ*sizeof(double),next_Sij, &next_Sij);
            C_DGEMM('N','N',nso,n1occ,n1occ,1.0,C1op[0],n1occ,Sij1p[0],n1occ,0.0,R[2*a]->pointer()[0],n1occ);
            Dx[a] = linalg::doublet(L[2*a], R[2*a], false, true);
            // Symmetrize the pseudodensity
            Dx[a]->add(Dx[a]->transpose());
            Dx[a]->scale(0.5);
        }
        for (int a = 0; a < nA; a++) {
            psio_address next_Sij= psio_get_address(PSIO_ZERO,(A + a) * (size_t) n2occ * n2occ * sizeof(double));
            psio_->read(PSIF_HESS,Sij_2,(char*)Sij2p[0], static_cast<size_t> (n2occ)*n2occ*sizeof(double),next_Sij, &next_Sij);
            C_DGEMM('N','N',nso,n2occ,n2occ,1.0,C2op[0],n2occ,Sij2p[0],n2occ,0.0,R[2*a+1]->pointer()[0],n2occ);
        }

        jk->compute();
        if(functional_->needs_xc()) {
            potential_->compute_Vx(Dx, Vx);
        }


        for (int a = 0; a < nA; a++) {

            // Add the alpha J contribution to G
            C_DGEMM('N','N',nso,n1occ,nso,1.0,J[2*a]->pointer()[0],nso,C1op[0],n1occ,0.0,Tp[0],n1occ);
            C_DGEMM('T','N',nmo,n1occ,nso,-1.0,C1p[0],nmo,Tp[0],n1occ,0.0,Up[0],n1occ);
            T->zero();

            C_DGEMM('N','N',nso,n1occ,nso,1.0,J[2*a+1]->pointer()[0],nso,C1op[0],n1occ,0.0,Tp[0],n1occ);
            C_DGEMM('T','N',nmo,n1occ,nso,-1.0,C1p[0],nmo,Tp[0],n1occ,1.0,Up[0],n1occ);

            if(functional_->needs_xc()) {
                // Symmetrize the result, just to be safe
                C_DGEMM('N','N',nso,n1occ,nso, 0.5,Vx[2*a]->pointer()[0],nso,C1op[0],n1occ,0.0,Tp[0],n1occ);
                C_DGEMM('T','N',nso,n1occ,nso, 0.5,Vx[2*a]->pointer()[0],nso,C1op[0],n1occ,1.0,Tp[0],n1occ);
                C_DGEMM('T','N',nmo,n1occ,nso,-2.0,C1p[0],nmo,Tp[0],n1occ,1.0,Up[0],n1occ);
            }

            // Subtract the K term from G
            if( Kscale) {
                T->zero();
                C_DGEMM('N','N',nso,n1occ,nso,1.0,K[2*a]->pointer()[0],nso,C1op[0],n1occ,0.0,Tp[0],n1occ);
                C_DGEMM('T','N',nmo,n1occ,nso,Kscale,C1p[0],nmo,Tp[0],n1occ,1.0,Up[0],n1occ);
            }

            psio_address next_Gpi = psio_get_address(PSIO_ZERO,(A + a) * (size_t) nmo * n1occ * sizeof(double));
            psio_->write(PSIF_HESS,G2pi_str,(char*)Up[0], static_cast<size_t> (nmo)*n1occ*sizeof(double),next_Gpi,&next_Gpi);
        }
    }
}

void USCFDeriv::VXC_deriv(std::shared_ptr<Matrix> C, 
                          std::shared_ptr<Matrix> Cocc,
                          int nso, int nocc, int nvir, bool alpha)
{
    // => XC Gradient <= //
    double** Cp  = C->pointer();  
    double** Cop = Cocc->pointer();
    size_t nmo = nocc + nvir;
    int natom = molecule_->natom();

    auto T = std::make_shared<Matrix>("T",nso,nocc);
    double** Tp = T->pointer();
    auto U = std::make_shared<Matrix>("Tempai",nmo,nocc);
    double** Up = U->pointer();

    auto VXCpi_str = (alpha) ? "VXCpi^A_a" : "VXCpi^A_b";

    if (functional_->needs_xc()) {
        // Write some placeholder data to PSIO, to get the sizing right
        psio_address next_VXCpi = PSIO_ZERO;

        for (int A = 0; A < 3 * natom; A++)
            psio_->write(PSIF_HESS,VXCpi_str,(char*)Up[0], static_cast<size_t> (nmo)*nocc*sizeof(double),next_VXCpi,&next_VXCpi);

        // For now we just compute all 3N matrices in one go.  If this becomes to burdensome
        // in terms of memory we can reevaluate and implement a batching mechanism instead.
        auto Vxc_matrices = potential_->compute_fock_derivatives();
        for(int a =0; a < 3*natom; ++a){
            // Transform from SO basis to pi
            C_DGEMM('N','N',nso,nocc,nso,1.0,Vxc_matrices[a]->pointer()[0],nso,Cop[0],nocc,0.0,Tp[0],nocc);
            C_DGEMM('T','N',nmo,nocc,nso,1.0,Cp[0],nmo,Tp[0],nocc,0.0,Up[0],nocc);
            next_VXCpi = psio_get_address(PSIO_ZERO,a * (size_t) nmo * nocc * sizeof(double));

            psio_->write(PSIF_HESS,VXCpi_str,(char*)Up[0], static_cast<size_t> (nmo)*nocc*sizeof(double),next_VXCpi,&next_VXCpi);
        }
    }
}

void USCFDeriv::assemble_Fock(int nocc, int nvir, bool alpha)
{
    // => Fpi <= //

    size_t nmo = nocc + nvir;
    int natom = molecule_->natom();
    auto Tpi = std::make_shared<Matrix>("Tpi",nmo,nocc);
    auto Fpi = std::make_shared<Matrix>("Fpi",nmo,nocc);
    double** Tpip = Tpi->pointer();
    double** Fpip = Fpi->pointer();

    psio_address next_Tpi = PSIO_ZERO;
    psio_address next_Vpi = PSIO_ZERO;
    psio_address next_Jpi = PSIO_ZERO;
    psio_address next_Fpi = PSIO_ZERO;
    psio_address next_Epi = PSIO_ZERO;
    psio_address next_VXCpi = PSIO_ZERO;


    auto T_str = (alpha) ? "Tpi^A_a" : "Tpi^A_b";
    auto V_str = (alpha) ? "Vpi^A_a" : "Vpi^A_b";
    auto G_str = (alpha) ? "Gpi^A_a" : "Gpi^A_b";
    auto E_str = (alpha) ? "Epi^A_a" : "Epi^A_b";
    auto VXC_str = (alpha) ? "VXCpi^A_a" : "VXCpi^A_b";
    auto F_str = (alpha) ? "Fpi^A_a" : "Fpi^A_b";

    // Add all of the same-spin components
    // opposite spin already included in G
    for (int A = 0; A < 3*natom; A++) {
        psio_->read(PSIF_HESS,T_str,(char*)Fpip[0], static_cast<size_t> (nmo) * nocc * sizeof(double),next_Tpi,&next_Tpi);
        psio_->read(PSIF_HESS,V_str,(char*)Tpip[0], static_cast<size_t> (nmo) * nocc * sizeof(double),next_Vpi,&next_Vpi);
        Fpi->add(Tpi);
        if (basisset_->has_ECP()) {
            psio_->read(PSIF_HESS,E_str,(char*)Tpip[0], static_cast<size_t> (nmo) * nocc * sizeof(double),next_Epi,&next_Epi);
            Fpi->add(Tpi);
        }
        psio_->read(PSIF_HESS,G_str,(char*)Tpip[0], static_cast<size_t> (nmo) * nocc * sizeof(double),next_Jpi,&next_Jpi);
        Fpi->add(Tpi);
        if (functional_->needs_xc()) {
            psio_->read(PSIF_HESS,VXC_str,(char*)Tpip[0], static_cast<size_t> (nmo) * nocc * sizeof(double),next_VXCpi,&next_VXCpi);
            Fpi->add(Tpi);
        }

        psio_->write(PSIF_HESS,F_str,(char*)Fpip[0], static_cast<size_t> (nmo) * nocc * sizeof(double),next_Fpi,&next_Fpi);
    }
}

void USCFDeriv::assemble_B(std::shared_ptr<Vector> eocc, int nocc, int nvir, bool alpha)
{
    // => Bpi <= //
    auto Tai = std::make_shared<Matrix>("T",nvir,nocc);
    auto Bai = std::make_shared<Matrix>("B",nvir,nocc);
    double** Taip = Tai->pointer();
    double** Baip = Bai->pointer();
    size_t nmo = nocc + nvir;
    int natom = molecule_->natom();

    double* eop = eocc->pointer();

    psio_address next_Fpi = PSIO_ZERO;
    psio_address next_Spi = PSIO_ZERO;
    psio_address next_G2pi = PSIO_ZERO;
    psio_address next_Bai = PSIO_ZERO;

    auto F_str = (alpha) ? "Fpi^A_a" : "Fpi^A_b";
    auto S_str = (alpha) ? "Spi^A_a" : "Spi^A_b";
    auto G_str = (alpha) ? "G2pi^A_a" : "G2pi^A_b";
    auto B_str = (alpha) ? "Bai^A_a" : "Bai^A_b";

    for (int A = 0; A < 3*natom; A++) {
        next_Fpi = psio_get_address(PSIO_ZERO,sizeof(double)*(A * (size_t) nmo * nocc + nocc * nocc));
        psio_->read(PSIF_HESS,F_str,(char*)Baip[0], static_cast<size_t> (nvir) * nocc * sizeof(double),next_Fpi,&next_Fpi);
        next_Spi = psio_get_address(PSIO_ZERO,sizeof(double)*(A * (size_t) nmo * nocc + nocc * nocc));
        psio_->read(PSIF_HESS,S_str,(char*)Taip[0], static_cast<size_t> (nvir) * nocc * sizeof(double),next_Spi,&next_Spi);
        for (int i = 0; i < nocc; i++)
            C_DAXPY(nvir,-eop[i],&Taip[0][i],nocc,&Baip[0][i],nocc);
        next_G2pi = psio_get_address(PSIO_ZERO,sizeof(double)*(A * (size_t) nmo * nocc + nocc * nocc));
        psio_->read(PSIF_HESS,G_str,(char*)Taip[0], static_cast<size_t> (nvir) * nocc * sizeof(double),next_G2pi,&next_G2pi);
        Bai->add(Tai);
        Bai->scale(-1.0);
        psio_->write(PSIF_HESS,B_str,(char*)Baip[0], static_cast<size_t> (nvir) * nocc * sizeof(double),next_Bai,&next_Bai);
    }
}


void USCFDeriv::assemble_U(int nocc, int nvir, bool alpha)
{
    // => Upi <= //
    size_t nmo = nocc + nvir;
    int natom = molecule_->natom();
    auto Upi = std::make_shared<Matrix>("U",nmo,nocc);
    double** Upqp = Upi->pointer();

    psio_address next_Spi = PSIO_ZERO;
    psio_address next_Uai = PSIO_ZERO;
    psio_address next_Upi = PSIO_ZERO;

    auto S_str   = (alpha)? "Sij^A_a" : "Sij^A_b";
    auto Uai_str = (alpha)? "Uai^A_a" : "Uai^A_b";
    auto Upi_str = (alpha)? "Upi^A_a" : "Upi^A_b";

    for (int A = 0; A < 3*natom; A++) {
        psio_->read(PSIF_HESS,S_str,(char*)Upqp[0], static_cast<size_t> (nocc) * nocc * sizeof(double),next_Spi,&next_Spi);
        C_DSCAL(nocc * (size_t) nocc,-0.5, Upqp[0], 1);
        psio_->read(PSIF_HESS,Uai_str,(char*)Upqp[nocc], static_cast<size_t> (nvir) * nocc * sizeof(double),next_Uai,&next_Uai);
        psio_->write(PSIF_HESS,Upi_str,(char*)Upqp[0], static_cast<size_t> (nmo) * nocc * sizeof(double),next_Upi,&next_Upi);
    }
}

void USCFDeriv::assemble_Q(std::shared_ptr<JK> jk,
                           std::shared_ptr<Matrix> C1, 
                           std::shared_ptr<Matrix> C1occ,
                           std::shared_ptr<Matrix> C2, 
                           std::shared_ptr<Matrix> C2occ,
                           int nso, int n1occ, int n2occ, int n1vir, bool alpha)
{
    // => Qpi <= //
    size_t nmo = n1occ + n1vir;
    int natom = molecule_->natom();
    size_t mem = 0.9 * memory_ / 8L;
    size_t per_A = 3L * nso * nso + 1L * n1occ * nso;
    size_t max_A = (mem / 2L) / per_A;

    double** C1p  = C1->pointer();  
    double** C1op = C1occ->pointer();
    double** C2p = C2->pointer();
    double** C2op = C2occ->pointer();

    max_A = (max_A > 3 * natom ? 3 * natom : max_A);
    double Kscale = functional_->x_alpha();
    std::vector<std::shared_ptr<Matrix> >& L = jk->C_left();
    std::vector<std::shared_ptr<Matrix> >& R = jk->C_right();
    const std::vector<std::shared_ptr<Matrix> >& J = jk->J();
    const std::vector<std::shared_ptr<Matrix> >& K = jk->K();
    std::vector<SharedMatrix> Dx, Vx;
    L.clear();
    R.clear();
    for (int a = 0; a < max_A; a++) {
        L.push_back(C1occ);
        R.push_back(std::make_shared<Matrix>("R",nso,n1occ));
        Dx.push_back(std::make_shared<Matrix>("Dx", nso,nso));
        Vx.push_back(std::make_shared<Matrix>("Vx", nso,nso));
        L.push_back(C2occ);
        R.push_back(std::make_shared<Matrix>("R",nso,n2occ));
    }

    auto U1pi = std::make_shared<Matrix>("Upi",nmo,n1occ);
    double** U1pip = U1pi->pointer();
    auto U2pi = std::make_shared<Matrix>("Upi",nmo,n2occ);
    double** U2pip = U2pi->pointer();
    auto T = std::make_shared<Matrix>("T",nso,n1occ);
    double** Tp = T->pointer();
    auto U = std::make_shared<Matrix>("T",nmo,n1occ);
    double** Up = U->pointer();
    
    auto Ustr_1 = (alpha) ? "Upi^A_a" : "Upi^A_b"; 
    auto Ustr_2 = (alpha) ? "Upi^A_b" : "Upi^A_a"; 
    auto Qstr = (alpha) ? "Qpi^A_a" :  "Qpi^A_b";

    for (int A = 0; A < 3 * natom; A+=max_A) {
        int nA = max_A;
        if (A + max_A >= 3 * natom) {
            nA = 3 * natom - A;
            L.resize(2*nA);
            R.resize(2*nA);
        }
        for (int a = 0; a < nA; a++) {
            psio_address next_Upi = psio_get_address(PSIO_ZERO,(A + a) * (size_t) nmo * n1occ * sizeof(double));
            psio_->read(PSIF_HESS,Ustr_1,(char*)U1pip[0], static_cast<size_t> (nmo)*n1occ*sizeof(double),next_Upi,&next_Upi);
            C_DGEMM('N','N',nso,n1occ,nmo,1.0,C1p[0],nmo,U1pip[0],n1occ,0.0,R[2*a]->pointer()[0],n1occ);
            Dx[a] = linalg::doublet(L[2*a], R[2*a], false, true);
            // Symmetrize the pseudodensity
            Dx[a]->add(Dx[a]->transpose());
            Dx[a]->scale(0.5);
        }
        for (int a = 0; a < nA; a++) {
            psio_address next_Upi = psio_get_address(PSIO_ZERO,(A + a) * (size_t) nmo * n2occ * sizeof(double));
            psio_->read(PSIF_HESS,Ustr_2,(char*)U2pip[0], static_cast<size_t> (nmo)*n2occ*sizeof(double),next_Upi,&next_Upi);
            C_DGEMM('N','N',nso,n2occ,nmo,1.0,C2p[0],nmo,U2pip[0],n2occ,0.0,R[2*a+1]->pointer()[0],n2occ);
        }

        jk->compute();
        if(functional_->needs_xc()) {
            potential_->compute_Vx(Dx, Vx);
        }

        for (int a = 0; a < nA; a++) {
            //C_DGEMM('N','N',nso,n1occ,nso, 4.0,J[a]->pointer()[0],nso,C1op[0],n1occ,0.0,Tp[0],n1occ);
            C_DGEMM('N','N',nso,n1occ,nso, 2.0,J[2*a]->pointer()[0],nso,C1op[0],n1occ,0.0,Tp[0],n1occ);
            C_DGEMM('N','N',nso,n1occ,nso, 2.0,J[2*a+1]->pointer()[0],nso,C1op[0],n1occ,1.0,Tp[0],n1occ);
            if(Kscale) {
                C_DGEMM('N','N',nso,n1occ,nso,-Kscale,K[2*a]->pointer()[0],nso,C1op[0],n1occ,1.0,Tp[0],n1occ);
                C_DGEMM('T','N',nso,n1occ,nso,-Kscale,K[2*a]->pointer()[0],nso,C1op[0],n1occ,1.0,Tp[0],n1occ);
            }
            if(functional_->needs_xc()) {
                // Symmetrize the result, just to be safe
                C_DGEMM('N','N',nso,n1occ,nso, 2.0,Vx[2*a]->pointer()[0],nso,C1op[0],n1occ,1.0,Tp[0],n1occ);
                C_DGEMM('T','N',nso,n1occ,nso, 2.0,Vx[2*a]->pointer()[0],nso,C1op[0],n1occ,1.0,Tp[0],n1occ);
            }
            C_DGEMM('T','N',nmo,n1occ,nso,1.0,C1p[0],nmo,Tp[0],n1occ,0.0,Up[0],n1occ);
            psio_address next_Qpi = psio_get_address(PSIO_ZERO,(A + a) * (size_t) nmo * n1occ * sizeof(double));
            psio_->write(PSIF_HESS,Qstr,(char*)Up[0], static_cast<size_t> (nmo)*n1occ*sizeof(double),next_Qpi,&next_Qpi);
        }
    }
}

void USCFDeriv::dipole_derivatives(std::shared_ptr<Matrix> C1, 
                           std::shared_ptr<Matrix> C1occ,
                           std::shared_ptr<Matrix> C2, 
                           std::shared_ptr<Matrix> C2occ,
                           std::shared_ptr<Matrix> Dt,
                           int nso, int n1occ, int n2occ, int n1vir)
{
    // => Dipole derivatives (for IR intensities) <= //

    size_t nmo = n1occ + n1vir;
    int natom = molecule_->natom();

    MintsHelper mints(basisset_);
    auto ao_dipole = mints.ao_dipole();
    Matrix mu1_x("mu X", nmo_, n1occ);
    Matrix mu1_y("mu Y", nmo_, n1occ);
    Matrix mu1_z("mu Z", nmo_, n1occ);
    Matrix mu2_x("mu X", nmo_, n2occ);
    Matrix mu2_y("mu Y", nmo_, n2occ);
    Matrix mu2_z("mu Z", nmo_, n2occ);
    mu1_x.transform(C1, ao_dipole[0], C1occ);
    mu1_y.transform(C1, ao_dipole[1], C1occ);
    mu1_z.transform(C1, ao_dipole[2], C1occ);
    mu2_x.transform(C2, ao_dipole[0], C2occ);
    mu2_y.transform(C2, ao_dipole[1], C2occ);
    mu2_z.transform(C2, ao_dipole[2], C2occ);
    // Start by computing the skeleton derivatives
    auto dipole_gradient = mints.dipole_grad(Dt);
    dipole_gradient->add(DipoleInt::nuclear_gradient_contribution(molecule_));
    double **pdip_grad = dipole_gradient->pointer();

    auto U1pi = std::make_shared<Matrix>("Upi",nmo,n1occ);
    double** U1pip = U1pi->pointer();
    auto U2pi = std::make_shared<Matrix>("Upi",nmo,n2occ);
    double** U2pip = U2pi->pointer();

    psio_address next_Upi = PSIO_ZERO;

    // TODO: Check this factor of 2
    for (int A = 0; A < 3*natom; A++) {
        psio_->read(PSIF_HESS,"Upi^A_a",(char*)U1pip[0], static_cast<size_t> (nmo) * n1occ * sizeof(double),next_Upi,&next_Upi);
        pdip_grad[A][0] += 2*mu1_x.vector_dot(U1pi);
        pdip_grad[A][1] += 2*mu1_y.vector_dot(U1pi);
        pdip_grad[A][2] += 2*mu1_z.vector_dot(U1pi);
    }
    for (int A = 0; A < 3*natom; A++) {
        psio_->read(PSIF_HESS,"Upi^A_b",(char*)U2pip[0], static_cast<size_t> (nmo) * n2occ * sizeof(double),next_Upi,&next_Upi);
        pdip_grad[A][0] += 2*mu2_x.vector_dot(U2pi);
        pdip_grad[A][1] += 2*mu2_y.vector_dot(U2pi);
        pdip_grad[A][2] += 2*mu2_z.vector_dot(U2pi);
    }
    uhf_wfn_->set_array_variable("SCF DIPOLE GRADIENT", dipole_gradient);
    uhf_wfn_->set_array_variable("CURRENT DIPOLE GRADIENT", dipole_gradient);
}

}  // namespace scfgrad
}  // namespace psi
