/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2019 The Psi4 Developers.
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
#include "psi4/libmints/sieve.h"
#include "psi4/libmints/cdsalclist.h"
#include "psi4/libfock/v.h"
#include "psi4/libfock/jk.h"
#include "psi4/libfock/apps.h"
#include "psi4/libfunctional/superfunctional.h"
#include "psi4/psifiles.h"
#include "psi4/lib3index/dftensor.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/vector.h"
#include "psi4/libmints/integral.h"
#include "psi4/libmints/mintshelper.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/libscf_solver/rhf.h"

#include <algorithm>

#include <sstream>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace psi;

namespace psi {
namespace scfgrad {

std::shared_ptr<Matrix> RSCFDeriv::hessian_response()
{
    // => Control Parameters <= //

    std::shared_ptr<Vector> eps     = epsilon_a_subset("AO","ALL");
    std::shared_ptr<Vector> eps_occ = epsilon_a_subset("AO","OCC");
    std::shared_ptr<Vector> eps_vir = epsilon_a_subset("AO","VIR");
    std::shared_ptr<Matrix> C    = Ca_subset("AO","ALL");
    std::shared_ptr<Matrix> Cocc = Ca_subset("AO","OCC");
    std::shared_ptr<Matrix> Cvir = Ca_subset("AO","VIR");
    std::shared_ptr<Matrix> Dt = Da_subset("AO");
    double** Dap = Dt->pointer();

    double** Cp  = C->pointer();
    double** Cop = Cocc->pointer();
    double** Cvp = Cvir->pointer();
    double*  ep  = eps->pointer();
    double*  eop = eps_occ->pointer();
    double*  evp = eps_vir->pointer();

    // => Sizing <= //

    int natom = molecule_->natom();
    int nso   = basisset_->nbf();
    int nocc  = eps_occ->dimpi()[0];
    int nvir  = eps_vir->dimpi()[0];
    int nmo   = C->colspi()[0];

    // => Target <= //

    auto response = std::make_shared<Matrix>("RHF Response",3*natom,3*natom);

    // => Response Utility File <= //

    psio_->open(PSIF_HESS,PSIO_OPEN_NEW);

    // => Spi <= //
    {
        // Overlap derivatives
        std::shared_ptr<OneBodyAOInt> Sint(integral_->ao_overlap(1));
        const double* buffer = Sint->buffer();

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
        for (int A = 0; A < 3*natom; A++) {
            psio_->write(PSIF_HESS,"Smi^A",(char*)Smixp[0], static_cast<size_t> (nso) * nocc * sizeof(double),next_Smi,&next_Smi);
        }
        for (int A = 0; A < 3*natom; A++) {
            psio_->write(PSIF_HESS,"Sai^A",(char*)Saip[0], static_cast<size_t> (nvir) * nocc * sizeof(double),next_Sai,&next_Sai);
        }
        for (int A = 0; A < 3*natom; A++) {
            psio_->write(PSIF_HESS,"Sij^A",(char*)Sijp[0], static_cast<size_t> (nocc) * nocc * sizeof(double),next_Sij,&next_Sij);
        }
        for (int A = 0; A < 3*natom; A++) {
            psio_->write(PSIF_HESS,"Spi^A",(char*)Spip[0], static_cast<size_t> (nmo) * nocc * sizeof(double),next_Spi,&next_Spi);
        }
        next_Smi = PSIO_ZERO;
        next_Sai = PSIO_ZERO;
        next_Sij = PSIO_ZERO;
        next_Spi = PSIO_ZERO;

        for (int A = 0; A < natom; A++) {
            Smix->zero();
            Smiy->zero();
            Smiz->zero();
            for (int P = 0; P < basisset_->nshell(); P++) {
                for (int Q = 0; Q < basisset_->nshell(); Q++) {
                    int aP = basisset_->shell(P).ncenter();
                    int aQ = basisset_->shell(Q).ncenter();
                    if ((aP != A && aQ != A) || aP == aQ) continue;
                    Sint->compute_shell_deriv1(P,Q);
                    int nP = basisset_->shell(P).nfunction();
                    int nQ = basisset_->shell(Q).nfunction();
                    int oP = basisset_->shell(P).function_index();
                    int oQ = basisset_->shell(Q).function_index();
                    const double* buffer2;
                    if (aP == A) {
                        // Px
                        buffer2 = buffer + 0 * nP * nQ;
                        for (int p = 0; p < nP; p++) {
                            for (int q = 0; q < nQ; q++) {
                                C_DAXPY(nocc,(*buffer2),Cop[q + oQ],1,Smixp[p + oP],1);
                                C_DAXPY(nocc,(*buffer2++), Cop[p + oP],1,Smixp[q + oQ],1);
                            }
                        }
                        // Py
                        buffer2 = buffer + 1 * nP * nQ;
                        for (int p = 0; p < nP; p++) {
                            for (int q = 0; q < nQ; q++) {
                                C_DAXPY(nocc,(*buffer2),Cop[q + oQ],1,Smiyp[p + oP],1);
                                C_DAXPY(nocc,(*buffer2++), Cop[p + oP],1,Smiyp[q + oQ],1);
                            }
                        }
                        // Pz
                        buffer2 = buffer + 2 * nP * nQ;
                        for (int p = 0; p < nP; p++) {
                            for (int q = 0; q < nQ; q++) {
                                C_DAXPY(nocc,(*buffer2),Cop[q + oQ],1,Smizp[p + oP],1);
                                C_DAXPY(nocc,(*buffer2++), Cop[p + oP],1,Smizp[q + oQ],1);
                            }
                        }
                    }
                    if (aQ == A) {
                        // Qx
                        buffer2 = buffer + 3 * nP * nQ;
                        for (int p = 0; p < nP; p++) {
                            for (int q = 0; q < nQ; q++) {
                                C_DAXPY(nocc,(*buffer2),Cop[q + oQ],1,Smixp[p + oP],1);
                                C_DAXPY(nocc,(*buffer2++), Cop[p + oP],1,Smixp[q + oQ],1);
                            }
                        }
                        // Qy
                        buffer2 = buffer + 4 * nP * nQ;
                        for (int p = 0; p < nP; p++) {
                            for (int q = 0; q < nQ; q++) {
                                C_DAXPY(nocc,(*buffer2),Cop[q + oQ],1,Smiyp[p + oP],1);
                                C_DAXPY(nocc,(*buffer2++), Cop[p + oP],1,Smiyp[q + oQ],1);
                            }
                        }
                        // Qz
                        buffer2 = buffer + 5 * nP * nQ;
                        for (int p = 0; p < nP; p++) {
                            for (int q = 0; q < nQ; q++) {
                                C_DAXPY(nocc,(*buffer2),Cop[q + oQ],1,Smizp[p + oP],1);
                                C_DAXPY(nocc,(*buffer2++), Cop[p + oP],1,Smizp[q + oQ],1);
                            }
                        }
                    }
                }
            }
            // Smi_x
            psio_->write(PSIF_HESS,"Smi^A",(char*)Smixp[0], static_cast<size_t> (nso) * nocc * sizeof(double),next_Smi,&next_Smi);
            // Smi_y
            psio_->write(PSIF_HESS,"Smi^A",(char*)Smiyp[0], static_cast<size_t> (nso) * nocc * sizeof(double),next_Smi,&next_Smi);
            // Smi_z
            psio_->write(PSIF_HESS,"Smi^A",(char*)Smizp[0], static_cast<size_t> (nso) * nocc * sizeof(double),next_Smi,&next_Smi);

            // Sai_x
            C_DGEMM('T','N',nvir,nocc,nso,0.5,Cvp[0],nvir,Smixp[0],nocc,0.0,Saip[0],nocc);
            psio_->write(PSIF_HESS,"Sai^A",(char*)Saip[0], static_cast<size_t> (nvir) * nocc * sizeof(double),next_Sai,&next_Sai);
            // Sai_y
            C_DGEMM('T','N',nvir,nocc,nso,0.5,Cvp[0],nvir,Smiyp[0],nocc,0.0,Saip[0],nocc);
            psio_->write(PSIF_HESS,"Sai^A",(char*)Saip[0], static_cast<size_t> (nvir) * nocc * sizeof(double),next_Sai,&next_Sai);
            // Sai_z
            C_DGEMM('T','N',nvir,nocc,nso,0.5,Cvp[0],nvir,Smizp[0],nocc,0.0,Saip[0],nocc);
            psio_->write(PSIF_HESS,"Sai^A",(char*)Saip[0], static_cast<size_t> (nvir) * nocc * sizeof(double),next_Sai,&next_Sai);

            // Sij_x
            C_DGEMM('T','N',nocc,nocc,nso,0.5,Cop[0],nocc,Smixp[0],nocc,0.0,Sijp[0],nocc);
            psio_->write(PSIF_HESS,"Sij^A",(char*)Sijp[0], static_cast<size_t> (nocc) * nocc * sizeof(double),next_Sij,&next_Sij);
            // Sij_y
            C_DGEMM('T','N',nocc,nocc,nso,0.5,Cop[0],nocc,Smiyp[0],nocc,0.0,Sijp[0],nocc);
            psio_->write(PSIF_HESS,"Sij^A",(char*)Sijp[0], static_cast<size_t> (nocc) * nocc * sizeof(double),next_Sij,&next_Sij);
            // Sij_z
            C_DGEMM('T','N',nocc,nocc,nso,0.5,Cop[0],nocc,Smizp[0],nocc,0.0,Sijp[0],nocc);
            psio_->write(PSIF_HESS,"Sij^A",(char*)Sijp[0], static_cast<size_t> (nocc) * nocc * sizeof(double),next_Sij,&next_Sij);

            // Spi_x
            C_DGEMM('T','N',nmo,nocc,nso,0.5,Cp[0],nmo,Smixp[0],nocc,0.0,Spip[0],nocc);
            psio_->write(PSIF_HESS,"Spi^A",(char*)Spip[0], static_cast<size_t> (nmo) * nocc * sizeof(double),next_Spi,&next_Spi);
            // Spi_y
            C_DGEMM('T','N',nmo,nocc,nso,0.5,Cp[0],nmo,Smiyp[0],nocc,0.0,Spip[0],nocc);
            psio_->write(PSIF_HESS,"Spi^A",(char*)Spip[0], static_cast<size_t> (nmo) * nocc * sizeof(double),next_Spi,&next_Spi);
            // Spi_z
            C_DGEMM('T','N',nmo,nocc,nso,0.5,Cp[0],nmo,Smizp[0],nocc,0.0,Spip[0],nocc);
            psio_->write(PSIF_HESS,"Spi^A",(char*)Spip[0], static_cast<size_t> (nmo) * nocc * sizeof(double),next_Spi,&next_Spi);
        }
    }

    // => Tpi <= //
    {
        // Kinetic derivatives
        std::shared_ptr<OneBodyAOInt> Tint(integral_->ao_kinetic(1));
        const double* buffer = Tint->buffer();

        auto Tmix = std::make_shared<Matrix>("Tmix",nso,nocc);
        auto Tmiy = std::make_shared<Matrix>("Tmiy",nso,nocc);
        auto Tmiz = std::make_shared<Matrix>("Tmiz",nso,nocc);
        double** Tmixp = Tmix->pointer();
        double** Tmiyp = Tmiy->pointer();
        double** Tmizp = Tmiz->pointer();

        auto Tpi = std::make_shared<Matrix>("Tpi",nmo,nocc);
        double** Tpip = Tpi->pointer();
        psio_address next_Tpi = PSIO_ZERO;


        for (int A = 0; A < natom; A++) {
            Tmix->zero();
            Tmiy->zero();
            Tmiz->zero();
            for (int P = 0; P < basisset_->nshell(); P++) {
                for (int Q = 0; Q < basisset_->nshell(); Q++) {
                    int aP = basisset_->shell(P).ncenter();
                    int aQ = basisset_->shell(Q).ncenter();
                    if ((aP != A && aQ != A) || aP == aQ) continue;
                    Tint->compute_shell_deriv1(P,Q);
                    int nP = basisset_->shell(P).nfunction();
                    int nQ = basisset_->shell(Q).nfunction();
                    int oP = basisset_->shell(P).function_index();
                    int oQ = basisset_->shell(Q).function_index();
                    const double* buffer2;
                    if (aP == A) {
                        // Px
                        buffer2 = buffer + 0 * nP * nQ;
                        for (int p = 0; p < nP; p++) {
                            for (int q = 0; q < nQ; q++) {
                                C_DAXPY(nocc,(*buffer2),Cop[q + oQ],1,Tmixp[p + oP],1);
                                C_DAXPY(nocc,(*buffer2++),Cop[p + oP],1,Tmixp[q + oQ],1);
                            }
                        }
                        // Py
                        buffer2 = buffer + 1 * nP * nQ;
                        for (int p = 0; p < nP; p++) {
                            for (int q = 0; q < nQ; q++) {
                                C_DAXPY(nocc,(*buffer2),Cop[q + oQ],1,Tmiyp[p + oP],1);
                                C_DAXPY(nocc,(*buffer2++),Cop[p + oP],1,Tmiyp[q + oQ],1);
                            }
                        }
                        // Pz
                        buffer2 = buffer + 2 * nP * nQ;
                        for (int p = 0; p < nP; p++) {
                            for (int q = 0; q < nQ; q++) {
                                C_DAXPY(nocc,(*buffer2),Cop[q + oQ],1,Tmizp[p + oP],1);
                                C_DAXPY(nocc,(*buffer2++),Cop[p + oP],1,Tmizp[q + oQ],1);
                            }
                        }
                    }
                    if (aQ == A) {
                        // Qx
                        buffer2 = buffer + 3 * nP * nQ;
                        for (int p = 0; p < nP; p++) {
                            for (int q = 0; q < nQ; q++) {
                                C_DAXPY(nocc,(*buffer2),Cop[q + oQ],1,Tmixp[p + oP],1);
                                C_DAXPY(nocc,(*buffer2++),Cop[p + oP],1,Tmixp[q + oQ],1);
                            }
                        }
                        // Qy
                        buffer2 = buffer + 4 * nP * nQ;
                        for (int p = 0; p < nP; p++) {
                            for (int q = 0; q < nQ; q++) {
                                C_DAXPY(nocc,(*buffer2),Cop[q + oQ],1,Tmiyp[p + oP],1);
                                C_DAXPY(nocc,(*buffer2++),Cop[p + oP],1,Tmiyp[q + oQ],1);
                            }
                        }
                        // Qz
                        buffer2 = buffer + 5 * nP * nQ;
                        for (int p = 0; p < nP; p++) {
                            for (int q = 0; q < nQ; q++) {
                                C_DAXPY(nocc,(*buffer2),Cop[q + oQ],1,Tmizp[p + oP],1);
                                C_DAXPY(nocc,(*buffer2++),Cop[p + oP],1,Tmizp[q + oQ],1);
                            }
                        }
                    }
                }
            }

            // Tpi_x
            C_DGEMM('T','N',nmo,nocc,nso,0.5,Cp[0],nmo,Tmixp[0],nocc,0.0,Tpip[0],nocc);
            psio_->write(PSIF_HESS,"Tpi^A",(char*)Tpip[0], static_cast<size_t> (nmo) * nocc * sizeof(double),next_Tpi,&next_Tpi);
            // Tpi_y
            C_DGEMM('T','N',nmo,nocc,nso,0.5,Cp[0],nmo,Tmiyp[0],nocc,0.0,Tpip[0],nocc);
            psio_->write(PSIF_HESS,"Tpi^A",(char*)Tpip[0], static_cast<size_t> (nmo) * nocc * sizeof(double),next_Tpi,&next_Tpi);
            // Tpi_z
            C_DGEMM('T','N',nmo,nocc,nso,0.5,Cp[0],nmo,Tmizp[0],nocc,0.0,Tpip[0],nocc);
            psio_->write(PSIF_HESS,"Tpi^A",(char*)Tpip[0], static_cast<size_t> (nmo) * nocc * sizeof(double),next_Tpi,&next_Tpi);
        }
    }


    // => Vpi <= //
    {
        // Potential derivatives
        std::shared_ptr<OneBodyAOInt> Vint(integral_->ao_potential(1));
        const double* buffer = Vint->buffer();

        auto Vmix = std::make_shared<Matrix>("Vmix",nso,nocc);
        auto Vmiy = std::make_shared<Matrix>("Vmiy",nso,nocc);
        auto Vmiz = std::make_shared<Matrix>("Vmiz",nso,nocc);
        double** Vmixp = Vmix->pointer();
        double** Vmiyp = Vmiy->pointer();
        double** Vmizp = Vmiz->pointer();

        auto Vpi = std::make_shared<Matrix>("Vpi",nmo,nocc);
        double** Vpip = Vpi->pointer();
        psio_address next_Vpi = PSIO_ZERO;


        for (int A = 0; A < natom; A++) {
            Vmix->zero();
            Vmiy->zero();
            Vmiz->zero();

            for (int P = 0; P < basisset_->nshell(); P++) {
                for (int Q = 0; Q < basisset_->nshell(); Q++) {
                    int aP = basisset_->shell(P).ncenter();
                    int aQ = basisset_->shell(Q).ncenter();
                    Vint->compute_shell_deriv1(P,Q);
                    int nP = basisset_->shell(P).nfunction();
                    int nQ = basisset_->shell(Q).nfunction();
                    int oP = basisset_->shell(P).function_index();
                    int oQ = basisset_->shell(Q).function_index();
                    const double* buf_x = &buffer[3 * A * nP * nQ + 0 * nP * nQ];
                    const double* buf_y = &buffer[3 * A * nP * nQ + 1 * nP * nQ];
                    const double* buf_z = &buffer[3 * A * nP * nQ + 2 * nP * nQ];

                    // Ax
                    for (int p = 0; p < nP; p++) {
                        for (int q = 0; q < nQ; q++) {
                            C_DAXPY(nocc,(*buf_x++),Cop[q + oQ],1,Vmixp[p + oP],1);
                        }
                    }

                    // Ay
                    for (int p = 0; p < nP; p++) {
                        for (int q = 0; q < nQ; q++) {
                            C_DAXPY(nocc,(*buf_y++),Cop[q + oQ],1,Vmiyp[p + oP],1);
                        }
                    }

                    // Az
                    for (int p = 0; p < nP; p++) {
                        for (int q = 0; q < nQ; q++) {
                            C_DAXPY(nocc,(*buf_z++),Cop[q + oQ],1,Vmizp[p + oP],1);
                        }
                    }
                }
            }

            // Vpi_x
            C_DGEMM('T','N',nmo,nocc,nso,1.0,Cp[0],nmo,Vmixp[0],nocc,0.0,Vpip[0],nocc);
            psio_->write(PSIF_HESS,"Vpi^A",(char*)Vpip[0], static_cast<size_t> (nmo) * nocc * sizeof(double),next_Vpi,&next_Vpi);
            // Vpi_y
            C_DGEMM('T','N',nmo,nocc,nso,1.0,Cp[0],nmo,Vmiyp[0],nocc,0.0,Vpip[0],nocc);
            psio_->write(PSIF_HESS,"Vpi^A",(char*)Vpip[0], static_cast<size_t> (nmo) * nocc * sizeof(double),next_Vpi,&next_Vpi);
            // Vpi_z
            C_DGEMM('T','N',nmo,nocc,nso,1.0,Cp[0],nmo,Vmizp[0],nocc,0.0,Vpip[0],nocc);
            psio_->write(PSIF_HESS,"Vpi^A",(char*)Vpip[0], static_cast<size_t> (nmo) * nocc * sizeof(double),next_Vpi,&next_Vpi);
        }
    }
#if 0
    ACS
    if (functional) {
        hessians_["Coulomb"] = jk_hessians["Coulomb"];
        if (functional->is_x_hybrid()) {
            hessians_["Exchange"] = jk_hessians["Exchange"];
            hessians_["Exchange"]->scale(-functional->x_alpha());
        }
        if (functional->is_x_lrc()) {
            hessians_["Exchange,LR"] = jk_hessians["Exchange,LR"];
            hessians_["Exchange,LR"]->scale(-functional->x_beta());
        }
    } else {
        hessians_["Coulomb"] = jk_hessians["Coulomb"];
        hessians_["Exchange"] = jk_hessians["Exchange"];
        hessians_["Exchange"]->scale(-1.0);
    }
#endif
    // => Jpi/Kpi <= //
    {
        // Figure out DFT functional info
        double Kscale = functional_->x_alpha();
        if (functional_->is_x_lrc())
            throw PSIEXCEPTION("Hessians for LRC functionals are not implemented yet.");

        size_t memory = 0.9 * memory_ / 8L;
        size_t max_a = memory / (3L * nso * nso);
        max_a = (max_a > 3 * natom ? 3 * natom : max_a);

        int natom = basisset_->molecule()->natom();

        std::vector<SharedMatrix> dGmats;
        std::vector<double**> pdG(3*natom);
        std::vector<bool> pert_incore(3*natom);
        for(int a = 0; a < max_a; ++a)
            dGmats.push_back(std::make_shared<Matrix>("G derivative contribution", nso, nso));

        auto Gpi = std::make_shared<Matrix>("MO G Deriv", nmo, nocc);
        double**pGpi = Gpi->pointer();
        psio_address next_Gpi = PSIO_ZERO;
        // Write some junk for now, to set the sizing for PSIO
        for(int pert = 0; pert < 3*natom; ++pert){
            psio_->write(PSIF_HESS,"Gpi^A",(char*)pGpi[0], static_cast<size_t> (nmo) * nocc * sizeof(double),next_Gpi,&next_Gpi);
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
                        const double* buffer = Pmnint->buffer();

                        for (int p = oP; p < oP+nP; p++) {
                            for (int m = oM; m < oM+nM; m++) {
                                for (int n = oN; n < oN+nN; n++) {
                                    Amnp[p][m*nso+n] += (*buffer++);
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

                    }
                }
            }
            // d[A] = Minv[A][B] c[B]  (factor of 2, to account for RHF)
            C_DGEMV('n', np, np, 2.0, PQp[0], np, cp, 1, 0.0, dp, 1);

            // B[B][m,n] = Minv[A][B] (A|mn)
            C_DGEMM('n','n', np, nso*nso, np, 1.0, PQp[0], np, Amnp[0], nso*nso, 0.0, Bmnp[0], nso*nso);

            // T[p][m,n] = B[p][r,n] D[m,r]
#pragma omp parallel for
            for(int p = 0; p < np; ++p)
                C_DGEMM('t', 'n', nso, nso, nso, 1.0, Dap[0], nso, Bmnp[p], nso, 0.0, pTmn[p], nso);


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

                        size_t stride = static_cast<size_t> (Pncart) * Qncart;

                        if(!pert_incore[Pcenter] && !pert_incore[Qcenter])
                            continue;

                        PQint->compute_shell_deriv1(P,0,Q,0);
                        const double* buffer = PQint->buffer();

                        auto *ptr = const_cast<double*>(buffer);

                        if(pert_incore[Pcenter]){
                            // J terms
                            // Px
                            C_DGEMV('n', nP, nQ, 1.0, ptr+0*stride, nQ, &dp[oQ], 1, 0.0, pTempP[0], 1);
                            C_DGEMV('t', nP, nso*nso, -1.0, Bmnp[oP], nso*nso, pTempP[0], 1, 1.0, pdG[Px][0], 1);
                            // Py
                            C_DGEMV('n', nP, nQ, 1.0, ptr+1*stride, nQ, &dp[oQ], 1, 0.0, pTempP[0], 1);
                            C_DGEMV('t', nP, nso*nso, -1.0, Bmnp[oP], nso*nso, pTempP[0], 1, 1.0, pdG[Py][0], 1);
                            // Pz
                            C_DGEMV('n', nP, nQ, 1.0, ptr+2*stride, nQ, &dp[oQ], 1, 0.0, pTempP[0], 1);
                            C_DGEMV('t', nP, nso*nso, -1.0, Bmnp[oP], nso*nso, pTempP[0], 1, 1.0, pdG[Pz][0], 1);

                            if (Kscale) {
                                // K terms
                                // Px
                                C_DGEMM('n', 'n', nP, nso*nso, nQ, 1.0, ptr+0*stride, nQ, pTmn[oQ], nso*nso, 0.0, pTmpPmn[0], nso*nso);
                                for(int p = 0; p < nP; ++p)
                                    C_DGEMM('N', 'N', nso, nso, nso, Kscale, Bmnp[p+oP], nso, pTmpPmn[p], nso, 1.0, pdG[Px][0], nso);
                                // Py
                                C_DGEMM('n', 'n', nP, nso*nso, nQ, 1.0, ptr+1*stride, nQ, pTmn[oQ], nso*nso, 0.0, pTmpPmn[0], nso*nso);
                                for(int p = 0; p < nP; ++p)
                                    C_DGEMM('N', 'N', nso, nso, nso, Kscale, Bmnp[p+oP], nso, pTmpPmn[p], nso, 1.0, pdG[Py][0], nso);
                                // Pz
                                C_DGEMM('n', 'n', nP, nso*nso, nQ, 1.0, ptr+2*stride, nQ, pTmn[oQ], nso*nso, 0.0, pTmpPmn[0], nso*nso);
                                for(int p = 0; p < nP; ++p)
                                    C_DGEMM('N', 'N', nso, nso, nso, Kscale, Bmnp[p+oP], nso, pTmpPmn[p], nso, 1.0, pdG[Pz][0], nso);
                            }

                        }
                        if(pert_incore[Qcenter]){
                            // J terms
                            // Qx
                            C_DGEMV('n', nP, nQ, 1.0, ptr+3*stride, nQ, &dp[oQ], 1, 0.0, pTempP[0], 1);
                            C_DGEMV('t', nP, nso*nso, -1.0, Bmnp[oP], nso*nso, pTempP[0], 1, 1.0, pdG[Qx][0], 1);
                            // Qy
                            C_DGEMV('n', nP, nQ, 1.0, ptr+4*stride, nQ, &dp[oQ], 1, 0.0, pTempP[0], 1);
                            C_DGEMV('t', nP, nso*nso, -1.0, Bmnp[oP], nso*nso, pTempP[0], 1, 1.0, pdG[Qy][0], 1);
                            // Qz
                            C_DGEMV('n', nP, nQ, 1.0, ptr+5*stride, nQ, &dp[oQ], 1, 0.0, pTempP[0], 1);
                            C_DGEMV('t', nP, nso*nso, -1.0, Bmnp[oP], nso*nso, pTempP[0], 1, 1.0, pdG[Qz][0], 1);

                            if (Kscale) {
                                // K terms
                                // Qx
                                C_DGEMM('n', 'n', nP, nso*nso, nQ, 1.0, ptr+3*stride, nQ, pTmn[oQ], nso*nso, 0.0, pTmpPmn[0], nso*nso);
                                for(int p = 0; p < nP; ++p)
                                    C_DGEMM('N', 'N', nso, nso, nso, Kscale, Bmnp[p+oP], nso, pTmpPmn[p], nso, 1.0, pdG[Qx][0], nso);
                                // Qy
                                C_DGEMM('n', 'n', nP, nso*nso, nQ, 1.0, ptr+4*stride, nQ, pTmn[oQ], nso*nso, 0.0, pTmpPmn[0], nso*nso);
                                for(int p = 0; p < nP; ++p)
                                    C_DGEMM('N', 'N', nso, nso, nso, Kscale, Bmnp[p+oP], nso, pTmpPmn[p], nso, 1.0, pdG[Qy][0], nso);
                                // Qz
                                C_DGEMM('n', 'n', nP, nso*nso, nQ, 1.0, ptr+5*stride, nQ, pTmn[oQ], nso*nso, 0.0, pTmpPmn[0], nso*nso);
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

                            size_t stride = static_cast<size_t> (Pncart) * Mncart * Nncart;

                            if(!pert_incore[Pcenter] &&
                                    !pert_incore[Mcenter] &&
                                    !pert_incore[Ncenter])
                                continue;

                            Pmnint->compute_shell_deriv1(P,0,M,N);
                            const double* buffer = Pmnint->buffer();

                            /*
                             * J terms have 2 contributions:
                             *      F^x[m][n] <- (P|mn)^x d[P]
                             * and
                             *      F^x[r][s] <- D[m][n] (P|mn)^x B[P][r,s]
                             * The second term factorizes into...
                             * ... Temp[P] = D[m][n] (P|mn)^x ...
                             * ... and then F^x[r][s] <- Temp[P] B[P][r,s]  (factor of 2 for RHF)
                             */

                            for(int x = 0; x < 9; ++ x){
                                size_t delta = 0L;
                                for(int p = 0; p < nP; ++p){
                                    double val = 0.0;
                                    for(int m = oM; m < nM+oM; ++m){
                                        for(int n = oN; n < nN+oN; ++n){
                                            val += Dap[m][n] * buffer[x*stride+delta];
                                            ++delta;
                                        }
                                    }
                                    pTempP[x][p] = val;
                                }
                            }

                            auto *ptr = const_cast<double*>(buffer);

                            if(pert_incore[Pcenter]){
                                // J Terms
                                size_t delta = 0L;
                                for(int p = oP; p < oP+nP; ++p){
                                    for(int m = oM; m < nM+oM; ++m){
                                        for(int n = oN; n < nN+oN; ++n){
                                            pdG[Px][m][n] += buffer[0*stride+delta] * dp[p];
                                            pdG[Py][m][n] += buffer[1*stride+delta] * dp[p];
                                            pdG[Pz][m][n] += buffer[2*stride+delta] * dp[p];
                                            ++delta;
                                        }
                                    }
                                }
                                C_DGEMV('t', nP, nso*nso, 2.0, Bmnp[oP], nso*nso, pTempP[0], 1, 1.0, pdG[Px][0], 1);
                                C_DGEMV('t', nP, nso*nso, 2.0, Bmnp[oP], nso*nso, pTempP[1], 1, 1.0, pdG[Py][0], 1);
                                C_DGEMV('t', nP, nso*nso, 2.0, Bmnp[oP], nso*nso, pTempP[2], 1, 1.0, pdG[Pz][0], 1);
                                // K Terms
                                if (Kscale) {
                                    for(int p = 0; p < nP; ++p)
                                        C_DGEMM('T', 'N', nN, nso, nM, -2 * Kscale, ptr+0*stride+p*nM*nN, nN, &pTmn[oP+p][oM*nso], nso, 1.0, pdG[Px][oN], nso);
                                    for(int p = 0; p < nP; ++p)
                                        C_DGEMM('T', 'N', nN, nso, nM, -2 * Kscale, ptr+1*stride+p*nM*nN, nN, &pTmn[oP+p][oM*nso], nso, 1.0, pdG[Py][oN], nso);
                                    for(int p = 0; p < nP; ++p)
                                        C_DGEMM('T', 'N', nN, nso, nM, -2 * Kscale, ptr+2*stride+p*nM*nN, nN, &pTmn[oP+p][oM*nso], nso, 1.0, pdG[Pz][oN], nso);
                                }
                            }
                            if(pert_incore[Mcenter]){
                                // J Terms
                                size_t delta = 0L;
                                for(int p = oP; p < oP+nP; ++p){
                                    for(int m = oM; m < nM+oM; ++m){
                                        for(int n = oN; n < nN+oN; ++n){
                                            pdG[mx][m][n] += buffer[3*stride+delta] * dp[p];
                                            pdG[my][m][n] += buffer[4*stride+delta] * dp[p];
                                            pdG[mz][m][n] += buffer[5*stride+delta] * dp[p];
                                            ++delta;
                                        }
                                    }
                                }
                                C_DGEMV('t', nP, nso*nso, 2.0, Bmnp[oP], nso*nso, pTempP[3], 1, 1.0, pdG[mx][0], 1);
                                C_DGEMV('t', nP, nso*nso, 2.0, Bmnp[oP], nso*nso, pTempP[4], 1, 1.0, pdG[my][0], 1);
                                C_DGEMV('t', nP, nso*nso, 2.0, Bmnp[oP], nso*nso, pTempP[5], 1, 1.0, pdG[mz][0], 1);
                                // K Terms
                                if (Kscale) {
                                    for(int p = 0; p < nP; ++p)
                                        C_DGEMM('T', 'N', nN, nso, nM, -2 * Kscale, ptr+3*stride+p*nM*nN, nN, &pTmn[oP+p][oM*nso], nso, 1.0, pdG[mx][oN], nso);
                                    for(int p = 0; p < nP; ++p)
                                        C_DGEMM('T', 'N', nN, nso, nM, -2 * Kscale, ptr+4*stride+p*nM*nN, nN, &pTmn[oP+p][oM*nso], nso, 1.0, pdG[my][oN], nso);
                                    for(int p = 0; p < nP; ++p)
                                        C_DGEMM('T', 'N', nN, nso, nM, -2 * Kscale, ptr+5*stride+p*nM*nN, nN, &pTmn[oP+p][oM*nso], nso, 1.0, pdG[mz][oN], nso);
                                }
                            }
                            if(pert_incore[Ncenter]){
                                // J Terms
                                size_t delta = 0L;
                                for(int p = oP; p < oP+nP; ++p){
                                    for(int m = oM; m < nM+oM; ++m){
                                        for(int n = oN; n < nN+oN; ++n){
                                            pdG[nx][m][n] += buffer[6*stride+delta] * dp[p];
                                            pdG[ny][m][n] += buffer[7*stride+delta] * dp[p];
                                            pdG[nz][m][n] += buffer[8*stride+delta] * dp[p];
                                            ++delta;
                                        }
                                    }
                                }
                                C_DGEMV('t', nP, nso*nso, 2.0, Bmnp[oP], nso*nso, pTempP[6], 1, 1.0, pdG[nx][0], 1);
                                C_DGEMV('t', nP, nso*nso, 2.0, Bmnp[oP], nso*nso, pTempP[7], 1, 1.0, pdG[ny][0], 1);
                                C_DGEMV('t', nP, nso*nso, 2.0, Bmnp[oP], nso*nso, pTempP[8], 1, 1.0, pdG[nz][0], 1);
                                // K Terms
                                if (Kscale) {
                                    for(int p = 0; p < nP; ++p)
                                        C_DGEMM('T', 'N', nN, nso, nM, -2 * Kscale, ptr+6*stride+p*nM*nN, nN, &pTmn[oP+p][oM*nso], nso, 1.0, pdG[nx][oN], nso);
                                    for(int p = 0; p < nP; ++p)
                                        C_DGEMM('T', 'N', nN, nso, nM, -2 * Kscale, ptr+7*stride+p*nM*nN, nN, &pTmn[oP+p][oM*nso], nso, 1.0, pdG[ny][oN], nso);
                                    for(int p = 0; p < nP; ++p)
                                        C_DGEMM('T', 'N', nN, nso, nM, -2 * Kscale, ptr+8*stride+p*nM*nN, nN, &pTmn[oP+p][oM*nso], nso, 1.0, pdG[nz][oN], nso);
                                }
                            }

                        }
                    }
                }

                for(int a = 0; a < nA; ++a){
                    // Symmetrize the derivative Fock contributions
                    SharedMatrix G = dGmats[a];
                    G->add(G->transpose());
                    Gpi->transform(C, G, Cocc);
                    Gpi->scale(0.5);
                    psio_->write(PSIF_HESS,"Gpi^A",(char*)pGpi[0], static_cast<size_t> (nmo) * nocc * sizeof(double),next_Gpi,&next_Gpi);
                }

            } // End loop over A batches

        }else{
            /*
             * The conventional integral algorithm
             */

            std::shared_ptr<TwoBodyAOInt> ints(integral_->eri(1));

            auto sieve = std::make_shared<ERISieve>(basisset_, 0.0);


            const std::vector<std::pair<int, int> >& shell_pairs = sieve->shell_pairs();
            size_t npairs = shell_pairs.size();
            size_t npairs2 = npairs * npairs;


            const double* buffer = ints->buffer();

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

                    if (!sieve->shell_significant(P,Q,R,S)) continue;

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

                    size_t stride = static_cast<size_t> (Pncart) * Qncart * Rncart * Sncart;

                    double Dpq, Drs, Dpr, Dqs, Dps, Dqr;
                    size_t delta;
                    double Ax, Ay, Az;
                    double Bx, By, Bz;
                    double Cx, Cy, Cz;
                    double Dx, Dy, Dz;

                    // => Coulomb Term <= //

                    Ax = 0.0; Ay = 0.0; Az = 0.0;
                    Bx = 0.0; By = 0.0; Bz = 0.0;
                    Cx = 0.0; Cy = 0.0; Cz = 0.0;
                    Dx = 0.0; Dy = 0.0; Dz = 0.0;
                    delta = 0L;
                    for (int p = Poff; p < Poff+Psize; p++) {
                        for (int q = Qoff; q < Qoff+Qsize; q++) {
                            for (int r = Roff; r < Roff+Rsize; r++) {
                                for (int s = Soff; s < Soff+Ssize; s++) {
                                    Dpq = Dap[p][q];
                                    Drs = Dap[r][s];
                                    Ax = prefactor * buffer[0 * stride + delta];
                                    Ay = prefactor * buffer[1 * stride + delta];
                                    Az = prefactor * buffer[2 * stride + delta];
                                    Cx = prefactor * buffer[3 * stride + delta];
                                    Cy = prefactor * buffer[4 * stride + delta];
                                    Cz = prefactor * buffer[5 * stride + delta];
                                    Dx = prefactor * buffer[6 * stride + delta];
                                    Dy = prefactor * buffer[7 * stride + delta];
                                    Dz = prefactor * buffer[8 * stride + delta];
                                    Bx = -(Ax + Cx + Dx);
                                    By = -(Ay + Cy + Dy);
                                    Bz = -(Az + Cz + Dz);

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
                                        Ax = prefactor * buffer[0 * stride + delta];
                                        Ay = prefactor * buffer[1 * stride + delta];
                                        Az = prefactor * buffer[2 * stride + delta];
                                        Cx = prefactor * buffer[3 * stride + delta];
                                        Cy = prefactor * buffer[4 * stride + delta];
                                        Cz = prefactor * buffer[5 * stride + delta];
                                        Dx = prefactor * buffer[6 * stride + delta];
                                        Dy = prefactor * buffer[7 * stride + delta];
                                        Dz = prefactor * buffer[8 * stride + delta];
                                        Bx = -(Ax + Cx + Dx);
                                        By = -(Ay + Cy + Dy);
                                        Bz = -(Az + Cz + Dz);

                                        Dpr = Dap[p][r];
                                        Dqs = Dap[q][s];
                                        Dps = Dap[p][s];
                                        Dqr = Dap[q][r];
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
                    Gpi->transform(C, G, Cocc);
                    Gpi->scale(0.5);
                    psio_->write(PSIF_HESS,"Gpi^A",(char*)pGpi[0], static_cast<size_t> (nmo) * nocc * sizeof(double),next_Gpi,&next_Gpi);
                }
            } // End loop over A batches

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
        if (functional_->is_x_lrc())
            throw PSIEXCEPTION("Hessians for LRC functionals are not implemented yet.");

        std::vector<std::shared_ptr<Matrix> >& L = jk->C_left();
        std::vector<std::shared_ptr<Matrix> >& R = jk->C_right();
        const std::vector<std::shared_ptr<Matrix> >& J = jk->J();
        const std::vector<std::shared_ptr<Matrix> >& K = jk->K();
        L.clear();
        R.clear();

        auto Sij = std::make_shared<Matrix>("Sij",nocc,nocc);
        double** Sijp = Sij->pointer();
        auto T = std::make_shared<Matrix>("T",nso,nocc);
        double** Tp = T->pointer();
        auto U = std::make_shared<Matrix>("Tempai",nmo,nocc);
        double** Up = U->pointer();

        // Write some placeholder data to PSIO, to get the sizing right
        psio_address next_Gpi = PSIO_ZERO;
        for (int A = 0; A < 3 * natom; A++)
            psio_->write(PSIF_HESS,"G2pi^A",(char*)Up[0], static_cast<size_t> (nmo)*nocc*sizeof(double),next_Gpi,&next_Gpi);

        std::vector<SharedMatrix> Dx, Vx;
        for (int A = 0; A < max_A; A++) {
            // Just pass C1 quantities in; this object doesn't respect symmetry anyway
            L.push_back(Cocc);
            R.push_back(std::make_shared<Matrix>("R",nso,nocc));
            Dx.push_back(std::make_shared<Matrix>("Dx", nso,nso));
            Vx.push_back(std::make_shared<Matrix>("Vx", nso,nso));
        }

        jk->print_header();

        for (int A = 0; A < 3 * natom; A+=max_A) {
            int nA = max_A;
            if (A + max_A >= 3 * natom) {
                nA = 3 * natom - A;
                L.resize(nA);
                R.resize(nA);
            }
            for (int a = 0; a < nA; a++) {
                psio_address next_Sij= psio_get_address(PSIO_ZERO,(A + a) * (size_t) nocc * nocc * sizeof(double));
                psio_->read(PSIF_HESS,"Sij^A",(char*)Sijp[0], static_cast<size_t> (nocc)*nocc*sizeof(double),next_Sij, &next_Sij);
                C_DGEMM('N','N',nso,nocc,nocc,1.0,Cop[0],nocc,Sijp[0],nocc,0.0,R[a]->pointer()[0],nocc);
                Dx[a] = linalg::doublet(L[a], R[a], false, true);
                // Symmetrize the pseudodensity
                Dx[a]->add(Dx[a]->transpose());
                Dx[a]->scale(0.5);
            }

            jk->compute();
            if(functional_->needs_xc()) {
                potential_->compute_Vx(Dx, Vx);
            }

            for (int a = 0; a < nA; a++) {
                // Add the 2J contribution to G
                C_DGEMM('N','N',nso,nocc,nso,1.0,J[a]->pointer()[0],nso,Cop[0],nocc,0.0,Tp[0],nocc);
                C_DGEMM('T','N',nmo,nocc,nso,-2.0,Cp[0],nmo,Tp[0],nocc,0.0,Up[0],nocc);

                if(functional_->needs_xc()) {
                    // Symmetrize the result, just to be safe
                    C_DGEMM('N','N',nso,nocc,nso, 0.5,Vx[a]->pointer()[0],nso,Cop[0],nocc,0.0,Tp[0],nocc);
                    C_DGEMM('T','N',nso,nocc,nso, 0.5,Vx[a]->pointer()[0],nso,Cop[0],nocc,1.0,Tp[0],nocc);
                    C_DGEMM('T','N',nmo,nocc,nso,-2.0,Cp[0],nmo,Tp[0],nocc,1.0,Up[0],nocc);
                }

                // Subtract the K term from G
                if( Kscale) {
                    C_DGEMM('N','N',nso,nocc,nso,1.0,K[a]->pointer()[0],nso,Cop[0],nocc,0.0,Tp[0],nocc);
                    C_DGEMM('T','N',nmo,nocc,nso,Kscale,Cp[0],nmo,Tp[0],nocc,1.0,Up[0],nocc);
                }

                psio_address next_Gpi = psio_get_address(PSIO_ZERO,(A + a) * (size_t) nmo * nocc * sizeof(double));
                psio_->write(PSIF_HESS,"G2pi^A",(char*)Up[0], static_cast<size_t> (nmo)*nocc*sizeof(double),next_Gpi,&next_Gpi);
            }
        }
    }

    // => XC Gradient <= //
    {
        auto T = std::make_shared<Matrix>("T",nso,nocc);
        double** Tp = T->pointer();
        auto U = std::make_shared<Matrix>("Tempai",nmo,nocc);
        double** Up = U->pointer();
        if (functional_->needs_xc()) {
            // Write some placeholder data to PSIO, to get the sizing right
            psio_address next_VXCpi = PSIO_ZERO;
            for (int A = 0; A < 3 * natom; A++)
                psio_->write(PSIF_HESS,"VXCpi^A",(char*)Up[0], static_cast<size_t> (nmo)*nocc*sizeof(double),next_VXCpi,&next_VXCpi);
            // For now we just compute all 3N matrices in one go.  If this becomes to burdensome
            // in terms of memory we can reevaluate and implement a batching mechanism instead.
            auto Vxc_matrices = potential_->compute_fock_derivatives();
            for(int a =0; a < 3*natom; ++a){
                // Transform from SO basis to pi
                C_DGEMM('N','N',nso,nocc,nso,1.0,Vxc_matrices[a]->pointer()[0],nso,Cop[0],nocc,0.0,Tp[0],nocc);
                C_DGEMM('T','N',nmo,nocc,nso,1.0,Cp[0],nmo,Tp[0],nocc,0.0,Up[0],nocc);
                next_VXCpi = psio_get_address(PSIO_ZERO,a * (size_t) nmo * nocc * sizeof(double));
                psio_->write(PSIF_HESS,"VXCpi^A",(char*)Up[0], static_cast<size_t> (nmo)*nocc*sizeof(double),next_VXCpi,&next_VXCpi);
            }
        }
    }

    // => Fpi <= //
    {
        auto Tpi = std::make_shared<Matrix>("Tpi",nmo,nocc);
        auto Fpi = std::make_shared<Matrix>("Fpi",nmo,nocc);
        double** Tpip = Tpi->pointer();
        double** Fpip = Fpi->pointer();

        psio_address next_Tpi = PSIO_ZERO;
        psio_address next_Vpi = PSIO_ZERO;
        psio_address next_Jpi = PSIO_ZERO;
        psio_address next_Fpi = PSIO_ZERO;
        psio_address next_VXCpi = PSIO_ZERO;

        for (int A = 0; A < 3*natom; A++) {
            psio_->read(PSIF_HESS,"Tpi^A",(char*)Fpip[0], static_cast<size_t> (nmo) * nocc * sizeof(double),next_Tpi,&next_Tpi);
            psio_->read(PSIF_HESS,"Vpi^A",(char*)Tpip[0], static_cast<size_t> (nmo) * nocc * sizeof(double),next_Vpi,&next_Vpi);
            Fpi->add(Tpi);
            psio_->read(PSIF_HESS,"Gpi^A",(char*)Tpip[0], static_cast<size_t> (nmo) * nocc * sizeof(double),next_Jpi,&next_Jpi);
            Fpi->add(Tpi);
            if (functional_->needs_xc()) {
                psio_->read(PSIF_HESS,"VXCpi^A",(char*)Tpip[0], static_cast<size_t> (nmo) * nocc * sizeof(double),next_VXCpi,&next_VXCpi);
                Fpi->add(Tpi);
            }
            psio_->write(PSIF_HESS,"Fpi^A",(char*)Fpip[0], static_cast<size_t> (nmo) * nocc * sizeof(double),next_Fpi,&next_Fpi);
        }
    }

    // => Bpi <= //
    {
        auto Tai = std::make_shared<Matrix>("T",nvir,nocc);
        auto Bai = std::make_shared<Matrix>("B",nvir,nocc);
        double** Taip = Tai->pointer();
        double** Baip = Bai->pointer();

        psio_address next_Fpi = PSIO_ZERO;
        psio_address next_Spi = PSIO_ZERO;
        psio_address next_G2pi = PSIO_ZERO;
        psio_address next_Bai = PSIO_ZERO;

        for (int A = 0; A < 3*natom; A++) {
            next_Fpi = psio_get_address(PSIO_ZERO,sizeof(double)*(A * (size_t) nmo * nocc + nocc * nocc));
            psio_->read(PSIF_HESS,"Fpi^A",(char*)Baip[0], static_cast<size_t> (nvir) * nocc * sizeof(double),next_Fpi,&next_Fpi);
            next_Spi = psio_get_address(PSIO_ZERO,sizeof(double)*(A * (size_t) nmo * nocc + nocc * nocc));
            psio_->read(PSIF_HESS,"Spi^A",(char*)Taip[0], static_cast<size_t> (nvir) * nocc * sizeof(double),next_Spi,&next_Spi);
            for (int i = 0; i < nocc; i++)
                C_DAXPY(nvir,-eop[i],&Taip[0][i],nocc,&Baip[0][i],nocc);
            next_G2pi = psio_get_address(PSIO_ZERO,sizeof(double)*(A * (size_t) nmo * nocc + nocc * nocc));
            psio_->read(PSIF_HESS,"G2pi^A",(char*)Taip[0], static_cast<size_t> (nvir) * nocc * sizeof(double),next_G2pi,&next_G2pi);
            Bai->add(Tai);
            Bai->scale(-1.0);
            psio_->write(PSIF_HESS,"Bai^A",(char*)Baip[0], static_cast<size_t> (nvir) * nocc * sizeof(double),next_Bai,&next_Bai);
        }
    }

    // => CPHF (Uai) <= //
    {
        rhf_wfn_->set_jk(jk);

        psio_address next_Bai = PSIO_ZERO;
        psio_address next_Uai = PSIO_ZERO;

        auto T = std::make_shared<Matrix>("T",nvir,nocc);
        double** Tp = T->pointer();

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
                auto B = std::make_shared<Matrix>(ss.str(),nocc,nvir);
                psio_->read(PSIF_HESS,"Bai^A",(char*)Tp[0], static_cast<size_t> (nvir) * nocc * sizeof(double),next_Bai,&next_Bai);
                double** Bp = B->pointer();
                for (int i = 0; i < nocc; i++) {
                    C_DCOPY(nvir,&Tp[0][i],nocc,Bp[i],1);
                }
                b_vecs.push_back(B);
            }

            auto u_matrices = rhf_wfn_->cphf_solve(b_vecs, options_.get_double("SOLVER_CONVERGENCE"),
                                                   options_.get_int("SOLVER_MAXITER"), print_);

            // Result in x
            for (int a = 0; a < nA; a++) {
                std::stringstream ss;
                ss << "Perturbation " << a + A;
                u_matrices[a]->scale(-1);
                double** Xp = u_matrices[a]->pointer();
                for (int i = 0; i < nocc; i++) {
                    C_DCOPY(nvir,Xp[i],1,&Tp[0][i],nocc);
                }
                psio_->write(PSIF_HESS,"Uai^A",(char*)Tp[0], static_cast<size_t> (nvir) * nocc * sizeof(double),next_Uai,&next_Uai);
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
    double **pdip_grad = dipole_gradient->pointer();

    // => Upi <= //
    {
        auto Upi = std::make_shared<Matrix>("U",nmo,nocc);
        double** Upqp = Upi->pointer();

        psio_address next_Spi = PSIO_ZERO;
        psio_address next_Uai = PSIO_ZERO;
        psio_address next_Upi = PSIO_ZERO;

        for (int A = 0; A < 3*natom; A++) {
            psio_->read(PSIF_HESS,"Sij^A",(char*)Upqp[0], static_cast<size_t> (nocc) * nocc * sizeof(double),next_Spi,&next_Spi);
            C_DSCAL(nocc * (size_t) nocc,-0.5, Upqp[0], 1);
            psio_->read(PSIF_HESS,"Uai^A",(char*)Upqp[nocc], static_cast<size_t> (nvir) * nocc * sizeof(double),next_Uai,&next_Uai);
            psio_->write(PSIF_HESS,"Upi^A",(char*)Upqp[0], static_cast<size_t> (nmo) * nocc * sizeof(double),next_Upi,&next_Upi);
            pdip_grad[A][0] += 4*mu_x.vector_dot(Upi);
            pdip_grad[A][1] += 4*mu_y.vector_dot(Upi);
            pdip_grad[A][2] += 4*mu_z.vector_dot(Upi);
        }
    }
    rhf_wfn_->set_array_variable("SCF DIPOLE GRADIENT", dipole_gradient);
    rhf_wfn_->set_array_variable("CURRENT DIPOLE GRADIENT", dipole_gradient);

    // => Qpi <= //
    {
        double Kscale = functional_->x_alpha();
        std::vector<std::shared_ptr<Matrix> >& L = jk->C_left();
        std::vector<std::shared_ptr<Matrix> >& R = jk->C_right();
        const std::vector<std::shared_ptr<Matrix> >& J = jk->J();
        const std::vector<std::shared_ptr<Matrix> >& K = jk->K();
        std::vector<SharedMatrix> Dx, Vx;
        L.clear();
        R.clear();
        for (int a = 0; a < max_A; a++) {
            L.push_back(Cocc);
            R.push_back(std::make_shared<Matrix>("R",nso,nocc));
            Dx.push_back(std::make_shared<Matrix>("Dx", nso,nso));
            Vx.push_back(std::make_shared<Matrix>("Vx", nso,nso));
        }

        auto Upi = std::make_shared<Matrix>("Upi",nmo,nocc);
        double** Upip = Upi->pointer();
        auto T = std::make_shared<Matrix>("T",nso,nocc);
        double** Tp = T->pointer();
        auto U = std::make_shared<Matrix>("T",nmo,nocc);
        double** Up = U->pointer();

        for (int A = 0; A < 3 * natom; A+=max_A) {
            int nA = max_A;
            if (A + max_A >= 3 * natom) {
                nA = 3 * natom - A;
                L.resize(nA);
                R.resize(nA);
            }
            for (int a = 0; a < nA; a++) {
                psio_address next_Upi = psio_get_address(PSIO_ZERO,(A + a) * (size_t) nmo * nocc * sizeof(double));
                psio_->read(PSIF_HESS,"Upi^A",(char*)Upip[0], static_cast<size_t> (nmo)*nocc*sizeof(double),next_Upi,&next_Upi);
                C_DGEMM('N','N',nso,nocc,nmo,1.0,Cp[0],nmo,Upip[0],nocc,0.0,R[a]->pointer()[0],nocc);
                Dx[a] = linalg::doublet(L[a], R[a], false, true);
                // Symmetrize the pseudodensity
                Dx[a]->add(Dx[a]->transpose());
                Dx[a]->scale(0.5);
            }

            jk->compute();
            if(functional_->needs_xc()) {
                potential_->compute_Vx(Dx, Vx);
            }

            for (int a = 0; a < nA; a++) {
                C_DGEMM('N','N',nso,nocc,nso, 4.0,J[a]->pointer()[0],nso,Cop[0],nocc,0.0,Tp[0],nocc);
                if(Kscale) {
                    C_DGEMM('N','N',nso,nocc,nso,-Kscale,K[a]->pointer()[0],nso,Cop[0],nocc,1.0,Tp[0],nocc);
                    C_DGEMM('T','N',nso,nocc,nso,-Kscale,K[a]->pointer()[0],nso,Cop[0],nocc,1.0,Tp[0],nocc);
                }
                if(functional_->needs_xc()) {
                    // Symmetrize the result, just to be safe
                    C_DGEMM('N','N',nso,nocc,nso, 2.0,Vx[a]->pointer()[0],nso,Cop[0],nocc,1.0,Tp[0],nocc);
                    C_DGEMM('T','N',nso,nocc,nso, 2.0,Vx[a]->pointer()[0],nso,Cop[0],nocc,1.0,Tp[0],nocc);
                }
                C_DGEMM('T','N',nmo,nocc,nso,1.0,Cp[0],nmo,Tp[0],nocc,0.0,Up[0],nocc);
                psio_address next_Qpi = psio_get_address(PSIO_ZERO,(A + a) * (size_t) nmo * nocc * sizeof(double));
                psio_->write(PSIF_HESS,"Qpi^A",(char*)Up[0], static_cast<size_t> (nmo)*nocc*sizeof(double),next_Qpi,&next_Qpi);
            }
        }
    }
    jk.reset();

    // => Zipper <= //
    {
        size_t memory = 0.9 * memory_ / 8L;
        size_t npi = nmo * (size_t) nocc;
        size_t max_a = memory / (3L * npi);
        max_a = (max_a > 3 * natom ? 3 * natom : max_a);

        auto L = std::make_shared<Matrix>("L",max_a * nmo, nocc);
        auto R = std::make_shared<Matrix>("R",max_a * nmo, nocc);
        auto T = std::make_shared<Matrix>("T",max_a * nmo, nocc);
        double** Lp = L->pointer();
        double** Rp = R->pointer();
        double** Tp = T->pointer();

        double** Hp = response->pointer();

        // U^A F^B
        for (int A = 0; A < 3 * natom; A+=max_a) {
            int nA = (A + max_a >= 3 * natom ? 3 * natom - A : max_a);
            psio_address nextA = psio_get_address(PSIO_ZERO, A * npi * sizeof(double));
            psio_->read(PSIF_HESS,"Upi^A",(char*)Lp[0],sizeof(double) * nA * npi,nextA,&nextA);
            for (int B = 0; B < 3 * natom; B+=max_a) {
                int nB = (B + max_a >= 3 * natom ? 3 * natom - B : max_a);
                psio_address nextB = psio_get_address(PSIO_ZERO, B * npi * sizeof(double));
                psio_->read(PSIF_HESS,"Fpi^A",(char*)Rp[0],sizeof(double) * nB * npi,nextB,&nextB);
                for (int a = 0; a < nA; a++) {
                    for (int b = 0; b < nB; b++) {
                        Hp[A + a][B + b] += 4.0 * C_DDOT(npi,Lp[0] + a * npi,1,Rp[0] + b * npi,1);
                    }
                }
            }
        }

        // F^A U^B
        for (int A = 0; A < 3 * natom; A+=max_a) {
            int nA = (A + max_a >= 3 * natom ? 3 * natom - A : max_a);
            psio_address nextA = psio_get_address(PSIO_ZERO, A * npi * sizeof(double));
            psio_->read(PSIF_HESS,"Fpi^A",(char*)Lp[0],sizeof(double) * nA * npi,nextA,&nextA);
            for (int B = 0; B < 3 * natom; B+=max_a) {
                int nB = (B + max_a >= 3 * natom ? 3 * natom - B : max_a);
                psio_address nextB = psio_get_address(PSIO_ZERO, B * npi * sizeof(double));
                psio_->read(PSIF_HESS,"Upi^A",(char*)Rp[0],sizeof(double) * nB * npi,nextB,&nextB);
                for (int a = 0; a < nA; a++) {
                    for (int b = 0; b < nB; b++) {
                        Hp[A + a][B + b] += 4.0 * C_DDOT(npi,Lp[0] + a * npi,1,Rp[0] + b * npi,1);
                    }
                }

            }
        }


        // U^A U^B
        // N.B. We use the relationship U^a_ia = -U^a_ai - S^a_ai
        psio_address junk;
        for (int A = 0; A < 3 * natom; A+=max_a) {
            int nA = (A + max_a >= 3 * natom ? 3 * natom - A : max_a);
            psio_address nextA = psio_get_address(PSIO_ZERO, A * npi * sizeof(double));
            psio_->read(PSIF_HESS,"Upi^A",(char*)Lp[0],sizeof(double) * nA * npi,nextA,&junk);
            psio_->read(PSIF_HESS,"Spi^A",(char*)Tp[0],sizeof(double) * nA * npi,nextA,&junk);
            L->add(T);
            for (int i = 0; i < nocc; i++)
                C_DSCAL(static_cast<size_t> (nA)*nmo,eop[i],&Lp[0][i],nocc);
            for (int B = 0; B < 3 * natom; B+=max_a) {
                int nB = (B + max_a >= 3 * natom ? 3 * natom - B : max_a);
                psio_address nextB = psio_get_address(PSIO_ZERO, B * npi * sizeof(double));
                psio_->read(PSIF_HESS,"Upi^A",(char*)Rp[0],sizeof(double) * nB * npi,nextB,&junk);
                psio_->read(PSIF_HESS,"Spi^A",(char*)Tp[0],sizeof(double) * nB * npi,nextB,&junk);
                R->add(T);
                for (int a = 0; a < nA; a++) {
                    for (int b = 0; b < nB; b++) {
                        Hp[A + a][B + b] -= 2.0 * C_DDOT(npi,Lp[0] + a * npi,1,Rp[0] + b * npi,1);
                        Hp[B + b][A + a] -= 2.0 * C_DDOT(npi,Lp[0] + a * npi,1,Rp[0] + b * npi,1);
                    }
                }
            }
        }

        // S^A S^B
        for (int A = 0; A < 3 * natom; A+=max_a) {
            int nA = (A + max_a >= 3 * natom ? 3 * natom - A : max_a);
            psio_address nextA = psio_get_address(PSIO_ZERO, A * npi * sizeof(double));
            psio_->read(PSIF_HESS,"Spi^A",(char*)Lp[0],sizeof(double) * nA * npi,nextA,&nextA);
            for (int i = 0; i < nocc; i++)
                C_DSCAL(static_cast<size_t> (nA)*nmo,eop[i],&Lp[0][i],nocc);
            for (int B = 0; B < 3 * natom; B+=max_a) {
                int nB = (B + max_a >= 3 * natom ? 3 * natom - B : max_a);
                psio_address nextB = psio_get_address(PSIO_ZERO, B * npi * sizeof(double));
                psio_->read(PSIF_HESS,"Spi^A",(char*)Rp[0],sizeof(double) * nB * npi,nextB,&nextB);
                for (int a = 0; a < nA; a++) {
                    for (int b = 0; b < nB; b++) {
                        Hp[A + a][B + b] += 2.0 * C_DDOT(npi,Lp[0] + a * npi,1,Rp[0] + b * npi,1);
                        Hp[B + b][A + a] += 2.0 * C_DDOT(npi,Lp[0] + a * npi,1,Rp[0] + b * npi,1);
                    }
                }
            }
        }

        // U^A U^B \epsilon
        for (int A = 0; A < 3 * natom; A+=max_a) {
            int nA = (A + max_a >= 3 * natom ? 3 * natom - A : max_a);
            psio_address nextA = psio_get_address(PSIO_ZERO, A * npi * sizeof(double));
            psio_->read(PSIF_HESS,"Upi^A",(char*)Lp[0],sizeof(double) * nA * npi,nextA,&nextA);
            double* Tp = Lp[0];
            for (int a = 0; a < nA; a++) {
                for (int p = 0; p < nmo; p++) {
                    C_DSCAL(nocc,ep[p],Tp,1);
                    Tp += nocc;
                }
            }
            for (int B = 0; B < 3 * natom; B+=max_a) {
                int nB = (B + max_a >= 3 * natom ? 3 * natom - B : max_a);
                psio_address nextB = psio_get_address(PSIO_ZERO, B * npi * sizeof(double));
                psio_->read(PSIF_HESS,"Upi^A",(char*)Rp[0],sizeof(double) * nB * npi,nextB,&nextB);
                for (int a = 0; a < nA; a++) {
                    for (int b = 0; b < nB; b++) {
                        Hp[A + a][B + b] += 4.0 * C_DDOT(npi,Lp[0] + a * npi,1,Rp[0] + b * npi,1);
                    }
                }
            }
        }

        // U^A Q^B
        for (int A = 0; A < 3 * natom; A+=max_a) {
            int nA = (A + max_a >= 3 * natom ? 3 * natom - A : max_a);
            psio_address nextA = psio_get_address(PSIO_ZERO, A * npi * sizeof(double));
            psio_->read(PSIF_HESS,"Upi^A",(char*)Lp[0],sizeof(double) * nA * npi,nextA,&nextA);
            for (int B = 0; B < 3 * natom; B+=max_a) {
                int nB = (B + max_a >= 3 * natom ? 3 * natom - B : max_a);
                psio_address nextB = psio_get_address(PSIO_ZERO, B * npi * sizeof(double));
                psio_->read(PSIF_HESS,"Qpi^A",(char*)Rp[0],sizeof(double) * nB * npi,nextB,&nextB);
                for (int a = 0; a < nA; a++) {
                    for (int b = 0; b < nB; b++) {
                        Hp[A + a][B + b] += 4.0 * C_DDOT(npi,Lp[0] + a * npi,1,Rp[0] + b * npi,1);
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


}} // Namespaces
