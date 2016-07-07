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

#include <libmints/mints.h>
#include <libqt/qt.h>
#include <libpsio/psio.hpp>
#include <libpsio/psio.h>
#include <psi4-dec.h>
#include <libfock/v.h>
#include <libfock/jk.h>
#include <libfock/apps.h>
#include <libfunctional/superfunctional.h>
#include <psifiles.h>
#include "libmints/sieve.h"
#include "libmints/view.h"
#include "scf_grad.h"
#include "jk_grad.h"


#include <algorithm>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>

#include <sstream>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;
using namespace psi;
using namespace boost;

namespace psi {
namespace scfgrad {

boost::shared_ptr<Matrix> SCFGrad::rhf_hessian_response()
{
    // => Control Parameters <= //

    boost::shared_ptr<Vector> eps     = epsilon_a_subset("AO","ALL");
    boost::shared_ptr<Vector> eps_occ = epsilon_a_subset("AO","OCC");
    boost::shared_ptr<Vector> eps_vir = epsilon_a_subset("AO","VIR");
    boost::shared_ptr<Matrix> C    = Ca_subset("AO","ALL");
    boost::shared_ptr<Matrix> Cocc = Ca_subset("AO","OCC");
    boost::shared_ptr<Matrix> Cvir = Ca_subset("AO","VIR");
    boost::shared_ptr<Matrix> Dt = Da_subset("AO");
    double** Dtp = Dt->pointer();

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
    
    boost::shared_ptr<Matrix> response(new Matrix("RHF Response",3*natom,3*natom));

    // => Response Utility File <= //

    psio_->open(PSIF_HESS,PSIO_OPEN_NEW);

    // => Spi <= //
    {
        // Overlap derivatives
        boost::shared_ptr<OneBodyAOInt> Sint(integral_->ao_overlap(1));
        const double* buffer = Sint->buffer();   

        boost::shared_ptr<Matrix> Smix(new Matrix("Smix",nso,nocc));
        boost::shared_ptr<Matrix> Smiy(new Matrix("Smiy",nso,nocc));
        boost::shared_ptr<Matrix> Smiz(new Matrix("Smiz",nso,nocc));
        double** Smixp = Smix->pointer();
        double** Smiyp = Smiy->pointer();
        double** Smizp = Smiz->pointer();

        boost::shared_ptr<Matrix> Sai(new Matrix("Sai",nvir,nocc));
        double** Saip = Sai->pointer();
        boost::shared_ptr<Matrix> Sij(new Matrix("Sij",nocc,nocc));
        double** Sijp = Sij->pointer();
        boost::shared_ptr<Matrix> Spi(new Matrix("Spi",nmo,nocc));
        double** Spip = Spi->pointer();

        psio_address next_Sai = PSIO_ZERO;
        psio_address next_Sij = PSIO_ZERO;
        psio_address next_Smi = PSIO_ZERO;
        psio_address next_Spi = PSIO_ZERO;
        for (int A = 0; A < 3*natom; A++) {
            psio_->write(PSIF_HESS,"Smi^A",(char*)Smixp[0],nso * nocc * sizeof(double),next_Smi,&next_Smi);
        }
        for (int A = 0; A < 3*natom; A++) {
            psio_->write(PSIF_HESS,"Sai^A",(char*)Saip[0],nvir * nocc * sizeof(double),next_Sai,&next_Sai);
        }
        for (int A = 0; A < 3*natom; A++) {
            psio_->write(PSIF_HESS,"Sij^A",(char*)Sijp[0],nocc * nocc * sizeof(double),next_Sij,&next_Sij);
        }
        for (int A = 0; A < 3*natom; A++) {
            psio_->write(PSIF_HESS,"Spi^A",(char*)Spip[0],nmo * nocc * sizeof(double),next_Spi,&next_Spi);
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
            psio_->write(PSIF_HESS,"Smi^A",(char*)Smixp[0],nso * nocc * sizeof(double),next_Smi,&next_Smi);
            // Smi_y
            psio_->write(PSIF_HESS,"Smi^A",(char*)Smiyp[0],nso * nocc * sizeof(double),next_Smi,&next_Smi);
            // Smi_z
            psio_->write(PSIF_HESS,"Smi^A",(char*)Smizp[0],nso * nocc * sizeof(double),next_Smi,&next_Smi);

            // Sai_x
            C_DGEMM('T','N',nvir,nocc,nso,0.5,Cvp[0],nvir,Smixp[0],nocc,0.0,Saip[0],nocc);
            psio_->write(PSIF_HESS,"Sai^A",(char*)Saip[0],nvir * nocc * sizeof(double),next_Sai,&next_Sai);
            // Sai_y
            C_DGEMM('T','N',nvir,nocc,nso,0.5,Cvp[0],nvir,Smiyp[0],nocc,0.0,Saip[0],nocc);
            psio_->write(PSIF_HESS,"Sai^A",(char*)Saip[0],nvir * nocc * sizeof(double),next_Sai,&next_Sai);
            // Sai_z
            C_DGEMM('T','N',nvir,nocc,nso,0.5,Cvp[0],nvir,Smizp[0],nocc,0.0,Saip[0],nocc);
            psio_->write(PSIF_HESS,"Sai^A",(char*)Saip[0],nvir * nocc * sizeof(double),next_Sai,&next_Sai);

            // Sij_x
            C_DGEMM('T','N',nocc,nocc,nso,0.5,Cop[0],nocc,Smixp[0],nocc,0.0,Sijp[0],nocc);
            psio_->write(PSIF_HESS,"Sij^A",(char*)Sijp[0],nocc * nocc * sizeof(double),next_Sij,&next_Sij);
            // Sij_y
            C_DGEMM('T','N',nocc,nocc,nso,0.5,Cop[0],nocc,Smiyp[0],nocc,0.0,Sijp[0],nocc);
            psio_->write(PSIF_HESS,"Sij^A",(char*)Sijp[0],nocc * nocc * sizeof(double),next_Sij,&next_Sij);
            // Sij_z
            C_DGEMM('T','N',nocc,nocc,nso,0.5,Cop[0],nocc,Smizp[0],nocc,0.0,Sijp[0],nocc);
            psio_->write(PSIF_HESS,"Sij^A",(char*)Sijp[0],nocc * nocc * sizeof(double),next_Sij,&next_Sij);

            // Spi_x
            C_DGEMM('T','N',nmo,nocc,nso,0.5,Cp[0],nmo,Smixp[0],nocc,0.0,Spip[0],nocc);
            psio_->write(PSIF_HESS,"Spi^A",(char*)Spip[0],nmo * nocc * sizeof(double),next_Spi,&next_Spi);
            // Spi_y
            C_DGEMM('T','N',nmo,nocc,nso,0.5,Cp[0],nmo,Smiyp[0],nocc,0.0,Spip[0],nocc);
            psio_->write(PSIF_HESS,"Spi^A",(char*)Spip[0],nmo * nocc * sizeof(double),next_Spi,&next_Spi);
            // Spi_z
            C_DGEMM('T','N',nmo,nocc,nso,0.5,Cp[0],nmo,Smizp[0],nocc,0.0,Spip[0],nocc);
            psio_->write(PSIF_HESS,"Spi^A",(char*)Spip[0],nmo * nocc * sizeof(double),next_Spi,&next_Spi);
        }
    }

    // => Tpi <= //
    {
        // Kinetic derivatives
        boost::shared_ptr<OneBodyAOInt> Tint(integral_->ao_kinetic(1));
        const double* buffer = Tint->buffer();   

        boost::shared_ptr<Matrix> Tmix(new Matrix("Tmix",nso,nocc));
        boost::shared_ptr<Matrix> Tmiy(new Matrix("Tmiy",nso,nocc));
        boost::shared_ptr<Matrix> Tmiz(new Matrix("Tmiz",nso,nocc));
        double** Tmixp = Tmix->pointer();
        double** Tmiyp = Tmiy->pointer();
        double** Tmizp = Tmiz->pointer();

        boost::shared_ptr<Matrix> Tpi(new Matrix("Tpi",nmo,nocc));
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
            psio_->write(PSIF_HESS,"Tpi^A",(char*)Tpip[0],nmo * nocc * sizeof(double),next_Tpi,&next_Tpi);
            // Tpi_y
            C_DGEMM('T','N',nmo,nocc,nso,0.5,Cp[0],nmo,Tmiyp[0],nocc,0.0,Tpip[0],nocc);
            psio_->write(PSIF_HESS,"Tpi^A",(char*)Tpip[0],nmo * nocc * sizeof(double),next_Tpi,&next_Tpi);
            // Tpi_z
            C_DGEMM('T','N',nmo,nocc,nso,0.5,Cp[0],nmo,Tmizp[0],nocc,0.0,Tpip[0],nocc);
            psio_->write(PSIF_HESS,"Tpi^A",(char*)Tpip[0],nmo * nocc * sizeof(double),next_Tpi,&next_Tpi);
        }
    }

    
    // => Vpi <= //
    {
        // Potential derivatives
        boost::shared_ptr<OneBodyAOInt> Vint(integral_->ao_potential(1));
        const double* buffer = Vint->buffer();   

        boost::shared_ptr<Matrix> Vmix(new Matrix("Vmix",nso,nocc));
        boost::shared_ptr<Matrix> Vmiy(new Matrix("Vmiy",nso,nocc));
        boost::shared_ptr<Matrix> Vmiz(new Matrix("Vmiz",nso,nocc));
        double** Vmixp = Vmix->pointer();
        double** Vmiyp = Vmiy->pointer();
        double** Vmizp = Vmiz->pointer();

        boost::shared_ptr<Matrix> Vpi(new Matrix("Vpi",nmo,nocc));
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
            psio_->write(PSIF_HESS,"Vpi^A",(char*)Vpip[0],nmo * nocc * sizeof(double),next_Vpi,&next_Vpi);
            // Vpi_y
            C_DGEMM('T','N',nmo,nocc,nso,1.0,Cp[0],nmo,Vmiyp[0],nocc,0.0,Vpip[0],nocc);
            psio_->write(PSIF_HESS,"Vpi^A",(char*)Vpip[0],nmo * nocc * sizeof(double),next_Vpi,&next_Vpi);
            // Vpi_z
            C_DGEMM('T','N',nmo,nocc,nso,1.0,Cp[0],nmo,Vmizp[0],nocc,0.0,Vpip[0],nocc);
            psio_->write(PSIF_HESS,"Vpi^A",(char*)Vpip[0],nmo * nocc * sizeof(double),next_Vpi,&next_Vpi);
        }
    }
    // => Jpi/Kpi <= //
    {

        size_t memory = 0.9 * memory_ / 8L;
        size_t max_a = memory / (3L * nso * nso);
        max_a = (max_a > 3 * natom ? 3 * natom : max_a);

        boost::shared_ptr<TwoBodyAOInt> ints(integral_->eri(1));

        int natom = basisset_->molecule()->natom();

        std::vector<SharedMatrix> dGmats;
        std::vector<double**> pdG(3*natom);
        std::vector<bool> pert_incore(3*natom);
        for(int a = 0; a < max_a; ++a)
            dGmats.push_back(SharedMatrix(new Matrix("G derivative contribution", nso, nso)));

        SharedMatrix Gpi = SharedMatrix(new Matrix("MO G Deriv", nmo, nocc));
        double**pGpi = Gpi->pointer();
        psio_address next_Gpi = PSIO_ZERO;
        // Write some junk for now, to set the sizing for PSIO
        for(int pert = 0; pert < 3*natom; ++pert){
            psio_->write(PSIF_HESS,"Gpi^A",(char*)pGpi[0],nmo * nocc * sizeof(double),next_Gpi,&next_Gpi);
        }
        next_Gpi = PSIO_ZERO;

        boost::shared_ptr<ERISieve> sieve(new ERISieve(basisset_, 0.0));


        const std::vector<std::pair<int, int> >& shell_pairs = sieve->shell_pairs();
        size_t npairs = shell_pairs.size();
        size_t npairs2 = npairs * npairs;


        const double* buffer = ints->buffer();

        for (int A = 0; A < 3 * natom; A+=max_a) {
            int nA = (A + max_a >= 3 * natom ? 3 * natom - A : max_a);

            // Keep track of which centers are loaded into memory, so we know when to skip
            std::fill(pert_incore.begin(), pert_incore.end(), false);
            std::fill(pdG.begin(), pdG.end(), (double**)NULL);
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
                    boost::swap(P,Q);
                    boost::swap(Pcenter,Qcenter);
                }
                if (am3 < am4){
                    boost::swap(R,S);
                    boost::swap(Rcenter,Scenter);
                }
                if( (am1 + am2) > (am3 + am4) ){
                    boost::swap(P,R);
                    boost::swap(Q,S);
                    boost::swap(Pcenter,Rcenter);
                    boost::swap(Qcenter,Scenter);
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

                size_t stride = Pncart * Qncart * Rncart * Sncart;

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
                                Dpq = Dtp[p][q];
                                Drs = Dtp[r][s];
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

                Ax = 0.0; Ay = 0.0; Az = 0.0;
                Bx = 0.0; By = 0.0; Bz = 0.0;
                Cx = 0.0; Cy = 0.0; Cz = 0.0;
                Dx = 0.0; Dy = 0.0; Dz = 0.0;
                delta = 0L;
                prefactor *= -0.25;
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

                                Dpr = Dtp[p][r];
                                Dqs = Dtp[q][s];
                                Dps = Dtp[p][s];
                                Dqr = Dtp[q][r];
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
            } // End shell loops

            for(int a = 0; a < nA; ++a){
                // Symmetrize the derivative Fock contributions
                SharedMatrix G = dGmats[a];
                G->add(G->transpose());
                Gpi->transform(C, G, Cocc);
                Gpi->scale(0.5);
                psio_->write(PSIF_HESS,"Gpi^A",(char*)pGpi[0],nmo * nocc * sizeof(double),next_Gpi,&next_Gpi);
            }
        } // End loop over A batches

    }

    boost::shared_ptr<JK> jk = JK::build_JK(basisset_, options_);
    jk->set_allow_desymmetrization(false);
    size_t mem = 0.9 * memory_ / 8L;
    size_t per_A = 3L * nso * nso + 1L * nocc * nso;
    size_t max_A = (mem / 2L) / per_A;
    max_A = (max_A > 3 * natom ? 3 * natom : max_A);
    jk->set_memory(mem);
    jk->initialize();

    // => J2pi/K2pi <= //
    {
        std::vector<boost::shared_ptr<Matrix> >& L = jk->C_left();  
        std::vector<boost::shared_ptr<Matrix> >& R = jk->C_right();  
        const std::vector<boost::shared_ptr<Matrix> >& J = jk->J();
        const std::vector<boost::shared_ptr<Matrix> >& K = jk->K();
        L.clear();
        R.clear();

        boost::shared_ptr<Matrix> Sij(new Matrix("Sij",nocc,nocc));
        double** Sijp = Sij->pointer();
        boost::shared_ptr<Matrix> T(new Matrix("T",nso,nocc));
        double** Tp = T->pointer();
        boost::shared_ptr<Matrix> U(new Matrix("Tempai",nmo,nocc));
        double** Up = U->pointer();

        // Write some placeholder data to PSIO, to get the sizing right
        psio_address next_Gpi = PSIO_ZERO;
        for (int A = 0; A < 3 * natom; A++)
            psio_->write(PSIF_HESS,"G2pi^A",(char*)Up[0],nmo*nocc*sizeof(double),next_Gpi,&next_Gpi);

        for (int A = 0; A < max_A; A++) {
            L.push_back(Cocc);
            R.push_back(boost::shared_ptr<Matrix>(new Matrix("R",nso,nocc)));
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
                psio_->read(PSIF_HESS,"Sij^A",(char*)Sijp[0],nocc*nocc*sizeof(double),next_Sij, &next_Sij);
                C_DGEMM('N','N',nso,nocc,nocc,1.0,Cop[0],nocc,Sijp[0],nocc,0.0,R[a]->pointer()[0],nocc);
            }

            jk->compute();
            for (int a = 0; a < nA; a++) {
                // Add the 2J contribution to G
                C_DGEMM('N','N',nso,nocc,nso,1.0,J[a]->pointer()[0],nso,Cop[0],nocc,0.0,Tp[0],nocc);
                C_DGEMM('T','N',nmo,nocc,nso,-2.0,Cp[0],nmo,Tp[0],nocc,0.0,Up[0],nocc);

                // Subtract the K term from G
                C_DGEMM('N','N',nso,nocc,nso,1.0,K[a]->pointer()[0],nso,Cop[0],nocc,0.0,Tp[0],nocc);
                C_DGEMM('T','N',nmo,nocc,nso,1.0,Cp[0],nmo,Tp[0],nocc,1.0,Up[0],nocc);
                psio_address next_Gpi = psio_get_address(PSIO_ZERO,(A + a) * (size_t) nmo * nocc * sizeof(double));
                psio_->write(PSIF_HESS,"G2pi^A",(char*)Up[0],nmo*nocc*sizeof(double),next_Gpi,&next_Gpi);
            }
        } 
    }

    // => Fpi <= //
    {
        boost::shared_ptr<Matrix> Tpi(new Matrix("Tpi",nmo,nocc));
        boost::shared_ptr<Matrix> Fpi(new Matrix("Fpi",nmo,nocc));
        double** Tpip = Tpi->pointer();
        double** Fpip = Fpi->pointer();

        psio_address next_Tpi = PSIO_ZERO;
        psio_address next_Vpi = PSIO_ZERO;
        psio_address next_Jpi = PSIO_ZERO;
        psio_address next_Fpi = PSIO_ZERO;

        for (int A = 0; A < 3*natom; A++) {
            psio_->read(PSIF_HESS,"Tpi^A",(char*)Fpip[0],nmo * nocc * sizeof(double),next_Tpi,&next_Tpi);
            psio_->read(PSIF_HESS,"Vpi^A",(char*)Tpip[0],nmo * nocc * sizeof(double),next_Vpi,&next_Vpi);
            Fpi->add(Tpi);
            psio_->read(PSIF_HESS,"Gpi^A",(char*)Tpip[0],nmo * nocc * sizeof(double),next_Jpi,&next_Jpi);
            Fpi->add(Tpi);
            psio_->write(PSIF_HESS,"Fpi^A",(char*)Fpip[0],nmo * nocc * sizeof(double),next_Fpi,&next_Fpi);
        }
    }

    // => Bpi <= //
    {
        boost::shared_ptr<Matrix> Tai(new Matrix("T",nvir,nocc));
        boost::shared_ptr<Matrix> Bai(new Matrix("B",nvir,nocc));
        double** Taip = Tai->pointer();
        double** Baip = Bai->pointer();

        psio_address next_Fpi = PSIO_ZERO;
        psio_address next_Spi = PSIO_ZERO;
        psio_address next_G2pi = PSIO_ZERO;
        psio_address next_Bai = PSIO_ZERO;

        for (int A = 0; A < 3*natom; A++) {
            next_Fpi = psio_get_address(PSIO_ZERO,sizeof(double)*(A * (size_t) nmo * nocc + nocc * nocc));
            psio_->read(PSIF_HESS,"Fpi^A",(char*)Baip[0],nvir * nocc * sizeof(double),next_Fpi,&next_Fpi);
            next_Spi = psio_get_address(PSIO_ZERO,sizeof(double)*(A * (size_t) nmo * nocc + nocc * nocc));
            psio_->read(PSIF_HESS,"Spi^A",(char*)Taip[0],nvir * nocc * sizeof(double),next_Spi,&next_Spi);
            for (int i = 0; i < nocc; i++)
                C_DAXPY(nvir,-eop[i],&Taip[0][i],nocc,&Baip[0][i],nocc);
            next_G2pi = psio_get_address(PSIO_ZERO,sizeof(double)*(A * (size_t) nmo * nocc + nocc * nocc));
            psio_->read(PSIF_HESS,"G2pi^A",(char*)Taip[0],nvir * nocc * sizeof(double),next_G2pi,&next_G2pi);
            Bai->add(Tai);
            Bai->scale(-1.0);
            psio_->write(PSIF_HESS,"Bai^A",(char*)Baip[0],nvir * nocc * sizeof(double),next_Bai,&next_Bai);
        }
    }    

    // => CPHF (Uai) <= //
    {
        SharedWavefunction wfn(new Wavefunction(options_));
        wfn->shallow_copy(this);
                
        boost::shared_ptr<RCPHF> cphf(new RCPHF(wfn, options_, false));
        cphf->set_jk(jk);

        std::map<std::string, SharedMatrix>& b = cphf->b();
        std::map<std::string, SharedMatrix>& x = cphf->x();

        psio_address next_Bai = PSIO_ZERO;    
        psio_address next_Uai = PSIO_ZERO;    

        boost::shared_ptr<Matrix> T(new Matrix("T",nvir,nocc));
        double** Tp = T->pointer();

        for (int A = 0; A < 3 * natom; A+=max_A) {
            int nA = max_A;
            if (A + max_A >= 3 * natom) {
                nA = 3 * natom - A;
            }
            
            x.clear();
            b.clear();

            // Fill b
            for (int a = 0; a < nA; a++) {
                psio_->read(PSIF_HESS,"Bai^A",(char*)Tp[0],nvir * nocc * sizeof(double),next_Bai,&next_Bai);
                std::stringstream ss;
                ss << "Perturbation " << a + A;
                boost::shared_ptr<Matrix> B(new Matrix(ss.str(),nocc,nvir));
                double** Bp = B->pointer();
                for (int i = 0; i < nocc; i++) {
                    C_DCOPY(nvir,&Tp[0][i],nocc,Bp[i],1);
                }
                b[ss.str()] = B;    
            }

            cphf->compute_energy();
            
            // Result in x
            for (int a = 0; a < nA; a++) {
                std::stringstream ss;
                ss << "Perturbation " << a + A;
                boost::shared_ptr<Matrix> X = x[ss.str()]; 
                double** Xp = X->pointer();
                for (int i = 0; i < nocc; i++) {
                    C_DCOPY(nvir,Xp[i],1,&Tp[0][i],nocc);
                }
                psio_->write(PSIF_HESS,"Uai^A",(char*)Tp[0],nvir * nocc * sizeof(double),next_Uai,&next_Uai);
            }
        } 
    }

    // => Upi <= //
    {
        boost::shared_ptr<Matrix> Upi(new Matrix("U",nmo,nocc));
        double** Upqp = Upi->pointer();

        psio_address next_Spi = PSIO_ZERO;    
        psio_address next_Uai = PSIO_ZERO;    
        psio_address next_Upi = PSIO_ZERO;    

        for (int A = 0; A < 3*natom; A++) {
            psio_->read(PSIF_HESS,"Sij^A",(char*)Upqp[0],nocc * nocc * sizeof(double),next_Spi,&next_Spi);
            C_DSCAL(nocc * (size_t) nocc,-0.5, Upqp[0], 1);
            psio_->read(PSIF_HESS,"Uai^A",(char*)Upqp[nocc],nvir * nocc * sizeof(double),next_Uai,&next_Uai);
            psio_->write(PSIF_HESS,"Upi^A",(char*)Upqp[0],nmo * nocc * sizeof(double),next_Upi,&next_Upi);
        }
    }

    // => Qpi <= //
    {
        std::vector<boost::shared_ptr<Matrix> >& L = jk->C_left();  
        std::vector<boost::shared_ptr<Matrix> >& R = jk->C_right();  
        const std::vector<boost::shared_ptr<Matrix> >& J = jk->J();
        const std::vector<boost::shared_ptr<Matrix> >& K = jk->K();
        L.clear();
        R.clear();
        for (int a = 0; a < max_A; a++) {
            L.push_back(Cocc);
            R.push_back(boost::shared_ptr<Matrix>(new Matrix("R",nso,nocc)));
        }

        boost::shared_ptr<Matrix> Upi(new Matrix("Upi",nmo,nocc));
        double** Upip = Upi->pointer();
        boost::shared_ptr<Matrix> T(new Matrix("T",nso,nocc));
        double** Tp = T->pointer();
        boost::shared_ptr<Matrix> U(new Matrix("T",nmo,nocc));
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
                psio_->read(PSIF_HESS,"Upi^A",(char*)Upip[0],nmo*nocc*sizeof(double),next_Upi,&next_Upi);
                C_DGEMM('N','N',nso,nocc,nmo,1.0,Cp[0],nmo,Upip[0],nocc,0.0,R[a]->pointer()[0],nocc);
            }

            jk->compute();
            for (int a = 0; a < nA; a++) {
                C_DGEMM('N','N',nso,nocc,nso, 4.0,J[a]->pointer()[0],nso,Cop[0],nocc,0.0,Tp[0],nocc);
                C_DGEMM('N','N',nso,nocc,nso,-1.0,K[a]->pointer()[0],nso,Cop[0],nocc,1.0,Tp[0],nocc);
                C_DGEMM('T','N',nso,nocc,nso,-1.0,K[a]->pointer()[0],nso,Cop[0],nocc,1.0,Tp[0],nocc);
                C_DGEMM('T','N',nmo,nocc,nso,1.0,Cp[0],nmo,Tp[0],nocc,0.0,Up[0],nocc);
                psio_address next_Qpi = psio_get_address(PSIO_ZERO,(A + a) * (size_t) nmo * nocc * sizeof(double));
                psio_->write(PSIF_HESS,"Qpi^A",(char*)Up[0],nmo*nocc*sizeof(double),next_Qpi,&next_Qpi);
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
        
        boost::shared_ptr<Matrix> L(new Matrix("L",max_a * nmo, nocc));
        boost::shared_ptr<Matrix> R(new Matrix("R",max_a * nmo, nocc));
        boost::shared_ptr<Matrix> T(new Matrix("T",max_a * nmo, nocc));
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
                C_DSCAL(nA*nmo,eop[i],&Lp[0][i],nocc);
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
                C_DSCAL(nA*nmo,eop[i],&Lp[0][i],nocc);
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
