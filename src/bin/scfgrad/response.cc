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

        boost::shared_ptr<Matrix> Spi(new Matrix("Spi",nmo,nocc));
        double** Spip = Spi->pointer();
        psio_address next_Spi = PSIO_ZERO;

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
                                C_DAXPY(nocc,(*buffer2++),Cop[q + oQ],1,Smixp[p + oP],1);
                            } 
                        } 
                        // Py
                        buffer2 = buffer + 1 * nP * nQ;                       
                        for (int p = 0; p < nP; p++) {
                            for (int q = 0; q < nQ; q++) {
                                C_DAXPY(nocc,(*buffer2++),Cop[q + oQ],1,Smiyp[p + oP],1);
                            } 
                        } 
                        // Pz
                        buffer2 = buffer + 2 * nP * nQ;                       
                        for (int p = 0; p < nP; p++) {
                            for (int q = 0; q < nQ; q++) {
                                C_DAXPY(nocc,(*buffer2++),Cop[q + oQ],1,Smizp[p + oP],1);
                            } 
                        } 
                    }
                    if (aQ == A) { 
                        // Qx
                        buffer2 = buffer + 3 * nP * nQ;                       
                        for (int p = 0; p < nP; p++) {
                            for (int q = 0; q < nQ; q++) {
                                C_DAXPY(nocc,(*buffer2++),Cop[q + oQ],1,Smixp[p + oP],1);
                            } 
                        } 
                        // Qy
                        buffer2 = buffer + 4 * nP * nQ;                       
                        for (int p = 0; p < nP; p++) {
                            for (int q = 0; q < nQ; q++) {
                                C_DAXPY(nocc,(*buffer2++),Cop[q + oQ],1,Smiyp[p + oP],1);
                            } 
                        } 
                        // Qz
                        buffer2 = buffer + 5 * nP * nQ;                       
                        for (int p = 0; p < nP; p++) {
                            for (int q = 0; q < nQ; q++) {
                                C_DAXPY(nocc,(*buffer2++),Cop[q + oQ],1,Smizp[p + oP],1);
                            } 
                        } 
                    }
                }
            }
            // Spi_x
            C_DGEMM('T','N',nmo,nocc,nso,1.0,Cp[0],nmo,Smixp[0],nocc,0.0,Spip[0],nocc);
            psio_->write(PSIF_HESS,"Spi^A",(char*)Spip[0],nmo * nocc * sizeof(double),next_Spi,&next_Spi);
            // Spi_y
            C_DGEMM('T','N',nmo,nocc,nso,1.0,Cp[0],nmo,Smiyp[0],nocc,0.0,Spip[0],nocc);
            psio_->write(PSIF_HESS,"Spi^A",(char*)Spip[0],nmo * nocc * sizeof(double),next_Spi,&next_Spi);
            // Spi_z
            C_DGEMM('T','N',nmo,nocc,nso,1.0,Cp[0],nmo,Smizp[0],nocc,0.0,Spip[0],nocc);
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
                                C_DAXPY(nocc,(*buffer2++),Cop[q + oQ],1,Tmixp[p + oP],1);
                            } 
                        } 
                        // Py
                        buffer2 = buffer + 1 * nP * nQ;                       
                        for (int p = 0; p < nP; p++) {
                            for (int q = 0; q < nQ; q++) {
                                C_DAXPY(nocc,(*buffer2++),Cop[q + oQ],1,Tmiyp[p + oP],1);
                            } 
                        } 
                        // Pz
                        buffer2 = buffer + 2 * nP * nQ;                       
                        for (int p = 0; p < nP; p++) {
                            for (int q = 0; q < nQ; q++) {
                                C_DAXPY(nocc,(*buffer2++),Cop[q + oQ],1,Tmizp[p + oP],1);
                            } 
                        } 
                    }
                    if (aQ == A) { 
                        // Qx
                        buffer2 = buffer + 3 * nP * nQ;                       
                        for (int p = 0; p < nP; p++) {
                            for (int q = 0; q < nQ; q++) {
                                C_DAXPY(nocc,(*buffer2++),Cop[q + oQ],1,Tmixp[p + oP],1);
                            } 
                        } 
                        // Qy
                        buffer2 = buffer + 4 * nP * nQ;                       
                        for (int p = 0; p < nP; p++) {
                            for (int q = 0; q < nQ; q++) {
                                C_DAXPY(nocc,(*buffer2++),Cop[q + oQ],1,Tmiyp[p + oP],1);
                            } 
                        } 
                        // Qz
                        buffer2 = buffer + 5 * nP * nQ;                       
                        for (int p = 0; p < nP; p++) {
                            for (int q = 0; q < nQ; q++) {
                                C_DAXPY(nocc,(*buffer2++),Cop[q + oQ],1,Tmizp[p + oP],1);
                            } 
                        } 
                    }
                }
            }
            // Tpi_x
            C_DGEMM('T','N',nmo,nocc,nso,1.0,Cp[0],nmo,Tmixp[0],nocc,0.0,Tpip[0],nocc);
            psio_->write(PSIF_HESS,"Tpi^A",(char*)Tpip[0],nmo * nocc * sizeof(double),next_Tpi,&next_Tpi);
            // Tpi_y
            C_DGEMM('T','N',nmo,nocc,nso,1.0,Cp[0],nmo,Tmiyp[0],nocc,0.0,Tpip[0],nocc);
            psio_->write(PSIF_HESS,"Tpi^A",(char*)Tpip[0],nmo * nocc * sizeof(double),next_Tpi,&next_Tpi);
            // Tpi_z
            C_DGEMM('T','N',nmo,nocc,nso,1.0,Cp[0],nmo,Tmizp[0],nocc,0.0,Tpip[0],nocc);
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
        // => Density matrix <= //

        boost::shared_ptr<Matrix> D(new Matrix("D",nso,nso));
        double** Dp = D->pointer();
        C_DGEMM('N','T',nso,nso,nocc,1.0,Cop[0],nocc,Cop[0],nocc,0.0,Dp[0],nso);

        boost::shared_ptr<TwoBodyAOInt> Iint(integral_->eri(1));
        const double* buffer = Iint->buffer();   

        boost::shared_ptr<Matrix> Jmix(new Matrix("Jmix",nso,nocc));
        boost::shared_ptr<Matrix> Jmiy(new Matrix("Jmiy",nso,nocc));
        boost::shared_ptr<Matrix> Jmiz(new Matrix("Jmiz",nso,nocc));
        double** Jmixp = Jmix->pointer();
        double** Jmiyp = Jmiy->pointer();
        double** Jmizp = Jmiz->pointer();

        boost::shared_ptr<Matrix> Kmix(new Matrix("Kmix",nso,nocc));
        boost::shared_ptr<Matrix> Kmiy(new Matrix("Kmiy",nso,nocc));
        boost::shared_ptr<Matrix> Kmiz(new Matrix("Kmiz",nso,nocc));
        double** Kmixp = Kmix->pointer();
        double** Kmiyp = Kmiy->pointer();
        double** Kmizp = Kmiz->pointer();

        boost::shared_ptr<Matrix> Tpi(new Matrix("Tpi",nmo,nocc));
        double** Tpip = Tpi->pointer();
        psio_address next_Jpi = PSIO_ZERO;
        psio_address next_Kpi = PSIO_ZERO;

        for (int A = 0; A < 3 * natom; A++) {
            psio_->write(PSIF_HESS,"Jpi^A",(char*)Tpip[0],nmo * nocc * sizeof(double),next_Jpi,&next_Jpi);
        }
        for (int A = 0; A < 3 * natom; A++) {
            psio_->write(PSIF_HESS,"Kpi^A",(char*)Tpip[0],nmo * nocc * sizeof(double),next_Kpi,&next_Kpi);
        }
        next_Jpi = PSIO_ZERO;
        next_Kpi = PSIO_ZERO;

        for (int A = 0; A < natom; A++) {
            Jmix->zero();
            Jmiy->zero();
            Jmiz->zero();
            Kmix->zero();
            Kmiy->zero();
            Kmiz->zero();

            for (int P = 0; P < basisset_->nshell(); P++) {
            for (int Q = 0; Q < basisset_->nshell(); Q++) {
            for (int R = 0; R < basisset_->nshell(); R++) {
            for (int S = 0; S < basisset_->nshell(); S++) {
                    int aP = basisset_->shell(P).ncenter();
                    int aQ = basisset_->shell(Q).ncenter();
                    int aR = basisset_->shell(R).ncenter();
                    int aS = basisset_->shell(S).ncenter();

                    if (!(aP == A | aQ == A | aR == A | aS == A)) continue;

                    Iint->compute_shell_deriv1(P,Q,R,S);

                    int nP = basisset_->shell(P).nfunction();
                    int nQ = basisset_->shell(Q).nfunction();
                    int nR = basisset_->shell(R).nfunction();
                    int nS = basisset_->shell(S).nfunction();

                    int cP = basisset_->shell(P).ncartesian();
                    int cQ = basisset_->shell(Q).ncartesian();
                    int cR = basisset_->shell(R).ncartesian();
                    int cS = basisset_->shell(S).ncartesian();

                    int oP = basisset_->shell(P).function_index();
                    int oQ = basisset_->shell(Q).function_index();
                    int oR = basisset_->shell(R).function_index();
                    int oS = basisset_->shell(S).function_index();

                    size_t cart_off = cP * cQ * cR * cS;

                    const double* Axp = buffer + 0 * cart_off;
                    const double* Ayp = buffer + 1 * cart_off;
                    const double* Azp = buffer + 2 * cart_off;
                    const double* Cxp = buffer + 3 * cart_off;
                    const double* Cyp = buffer + 4 * cart_off;
                    const double* Czp = buffer + 5 * cart_off;
                    const double* Dxp = buffer + 6 * cart_off;
                    const double* Dyp = buffer + 7 * cart_off;
                    const double* Dzp = buffer + 8 * cart_off;
            
                    // => J <= //

                    // A
                    if (aP == A) {
                        const double* A2xp = Axp;
                        const double* A2yp = Ayp;
                        const double* A2zp = Azp;

                        // Ax
                        for (int p = 0; p < nP; p++) {
                        for (int q = 0; q < nQ; q++) {
                            double Jval = 0.0;
                            for (int r = 0; r < nR; r++) {
                            for (int s = 0; s < nS; s++) {
                                Jval += Dp[r + oR][s + oS] * (*A2xp++); 
                            }}
                            C_DAXPY(nocc,Jval,Cop[q + oQ],1,Jmixp[p + oP],1);
                        }}

                        // Ay
                        for (int p = 0; p < nP; p++) {
                        for (int q = 0; q < nQ; q++) {
                            double Jval = 0.0;
                            for (int r = 0; r < nR; r++) {
                            for (int s = 0; s < nS; s++) {
                                Jval += Dp[r + oR][s + oS] * (*A2yp++); 
                            }}
                            C_DAXPY(nocc,Jval,Cop[q + oQ],1,Jmiyp[p + oP],1);
                        }}

                        // Az
                        for (int p = 0; p < nP; p++) {
                        for (int q = 0; q < nQ; q++) {
                            double Jval = 0.0;
                            for (int r = 0; r < nR; r++) {
                            for (int s = 0; s < nS; s++) {
                                Jval += Dp[r + oR][s + oS] * (*A2zp++); 
                            }}
                            C_DAXPY(nocc,Jval,Cop[q + oQ],1,Jmizp[p + oP],1);
                        }}
                    }

                    // B
                    if (aQ == A) {
                        const double* A2xp = Axp;
                        const double* A2yp = Ayp;
                        const double* A2zp = Azp;
                        const double* C2xp = Cxp;
                        const double* C2yp = Cyp;
                        const double* C2zp = Czp;
                        const double* D2xp = Dxp;
                        const double* D2yp = Dyp;
                        const double* D2zp = Dzp;

                        // Bx
                        for (int p = 0; p < nP; p++) {
                        for (int q = 0; q < nQ; q++) {
                            double Jval = 0.0;
                            for (int r = 0; r < nR; r++) {
                            for (int s = 0; s < nS; s++) {
                                Jval -= Dp[r + oR][s + oS] * ((*A2xp++) + (*C2xp++) + (*D2xp++)); 
                            }}
                            C_DAXPY(nocc,Jval,Cop[q + oQ],1,Jmixp[p + oP],1);
                        }}

                        // By
                        for (int p = 0; p < nP; p++) {
                        for (int q = 0; q < nQ; q++) {
                            double Jval = 0.0;
                            for (int r = 0; r < nR; r++) {
                            for (int s = 0; s < nS; s++) {
                                Jval -= Dp[r + oR][s + oS] * ((*A2yp++) + (*C2yp++) + (*D2yp++)); 
                            }}
                            C_DAXPY(nocc,Jval,Cop[q + oQ],1,Jmiyp[p + oP],1);
                        }}

                        // Bz
                        for (int p = 0; p < nP; p++) {
                        for (int q = 0; q < nQ; q++) {
                            double Jval = 0.0;
                            for (int r = 0; r < nR; r++) {
                            for (int s = 0; s < nS; s++) {
                                Jval -= Dp[r + oR][s + oS] * ((*A2zp++) + (*C2zp++) + (*D2zp++)); 
                            }}
                            C_DAXPY(nocc,Jval,Cop[q + oQ],1,Jmizp[p + oP],1);
                        }}
                    }

                    // C
                    if (aR == A) {
                        const double* C2xp = Cxp;
                        const double* C2yp = Cyp;
                        const double* C2zp = Czp;

                        // Cx
                        for (int p = 0; p < nP; p++) {
                        for (int q = 0; q < nQ; q++) {
                            double Jval = 0.0;
                            for (int r = 0; r < nR; r++) {
                            for (int s = 0; s < nS; s++) {
                                Jval += Dp[r + oR][s + oS] * (*C2xp++); 
                            }}
                            C_DAXPY(nocc,Jval,Cop[q + oQ],1,Jmixp[p + oP],1);
                        }}

                        // Cy
                        for (int p = 0; p < nP; p++) {
                        for (int q = 0; q < nQ; q++) {
                            double Jval = 0.0;
                            for (int r = 0; r < nR; r++) {
                            for (int s = 0; s < nS; s++) {
                                Jval += Dp[r + oR][s + oS] * (*C2yp++); 
                            }}
                            C_DAXPY(nocc,Jval,Cop[q + oQ],1,Jmiyp[p + oP],1);
                        }}

                        // Cz
                        for (int p = 0; p < nP; p++) {
                        for (int q = 0; q < nQ; q++) {
                            double Jval = 0.0;
                            for (int r = 0; r < nR; r++) {
                            for (int s = 0; s < nS; s++) {
                                Jval += Dp[r + oR][s + oS] * (*C2zp++); 
                            }}
                            C_DAXPY(nocc,Jval,Cop[q + oQ],1,Jmizp[p + oP],1);
                        }}
                    }

                    // D
                    if (aS == A) {
                        const double* D2xp = Dxp;
                        const double* D2yp = Dyp;
                        const double* D2zp = Dzp;

                        // Dx
                        for (int p = 0; p < nP; p++) {
                        for (int q = 0; q < nQ; q++) {
                            double Jval = 0.0;
                            for (int r = 0; r < nR; r++) {
                            for (int s = 0; s < nS; s++) {
                                Jval += Dp[r + oR][s + oS] * (*D2xp++); 
                            }}
                            C_DAXPY(nocc,Jval,Cop[q + oQ],1,Jmixp[p + oP],1);
                        }}

                        // Dy
                        for (int p = 0; p < nP; p++) {
                        for (int q = 0; q < nQ; q++) {
                            double Jval = 0.0;
                            for (int r = 0; r < nR; r++) {
                            for (int s = 0; s < nS; s++) {
                                Jval += Dp[r + oR][s + oS] * (*D2yp++); 
                            }}
                            C_DAXPY(nocc,Jval,Cop[q + oQ],1,Jmiyp[p + oP],1);
                        }}

                        // Dz
                        for (int p = 0; p < nP; p++) {
                        for (int q = 0; q < nQ; q++) {
                            double Jval = 0.0;
                            for (int r = 0; r < nR; r++) {
                            for (int s = 0; s < nS; s++) {
                                Jval += Dp[r + oR][s + oS] * (*D2zp++); 
                            }}
                            C_DAXPY(nocc,Jval,Cop[q + oQ],1,Jmizp[p + oP],1);
                        }}
                    }

                    // => K <= //

                    // A
                    if (aP == A) {

                        // Ax
                        for (int p = 0; p < nP; p++) {
                        for (int r = 0; r < nR; r++) {
                            double Kval = 0.0;
                            for (int q = 0; q < nQ; q++) {
                            const double* A2xp = Axp + p * nQ * nR * nS + q * nR * nS + r * nS;
                            for (int s = 0; s < nS; s++) {
                                Kval += Dp[p + oP][s + oS] * (*A2xp++); 
                            }}
                            C_DAXPY(nocc,Kval,Cop[r + oR],1,Kmixp[p + oP],1);
                        }}

                        // Ay
                        for (int p = 0; p < nP; p++) {
                        for (int r = 0; r < nR; r++) {
                            double Kval = 0.0;
                            for (int q = 0; q < nQ; q++) {
                            const double* A2yp = Ayp + p * nQ * nR * nS + q * nR * nS + r * nS;
                            for (int s = 0; s < nS; s++) {
                                Kval += Dp[p + oP][s + oS] * (*A2yp++); 
                            }}
                            C_DAXPY(nocc,Kval,Cop[r + oR],1,Kmiyp[p + oP],1);
                        }}

                        // Az
                        for (int p = 0; p < nP; p++) {
                        for (int r = 0; r < nR; r++) {
                            double Kval = 0.0;
                            for (int q = 0; q < nQ; q++) {
                            const double* A2zp = Azp + p * nQ * nR * nS + q * nR * nS + r * nS;
                            for (int s = 0; s < nS; s++) {
                                Kval += Dp[p + oP][s + oS] * (*A2zp++); 
                            }}
                            C_DAXPY(nocc,Kval,Cop[r + oR],1,Kmizp[p + oP],1);
                        }}
                    }

                    // B
                    if (aQ == A) {

                        // Bx
                        for (int p = 0; p < nP; p++) {
                        for (int r = 0; r < nR; r++) {
                            double Kval = 0.0;
                            for (int q = 0; q < nQ; q++) {
                            const double* A2xp = Axp + p * nQ * nR * nS + q * nR * nS + r * nS;
                            const double* C2xp = Cxp + p * nQ * nR * nS + q * nR * nS + r * nS;
                            const double* D2xp = Dxp + p * nQ * nR * nS + q * nR * nS + r * nS;
                            for (int s = 0; s < nS; s++) {
                                Kval -= Dp[p + oP][s + oS] * ((*A2xp++) + (*C2xp++) + (*D2xp++)); 
                            }}
                            C_DAXPY(nocc,Kval,Cop[r + oR],1,Kmixp[p + oP],1);
                        }}

                        // By
                        for (int p = 0; p < nP; p++) {
                        for (int r = 0; r < nR; r++) {
                            double Kval = 0.0;
                            for (int q = 0; q < nQ; q++) {
                            const double* A2yp = Ayp + p * nQ * nR * nS + q * nR * nS + r * nS;
                            const double* C2yp = Cyp + p * nQ * nR * nS + q * nR * nS + r * nS;
                            const double* D2yp = Dyp + p * nQ * nR * nS + q * nR * nS + r * nS;
                            for (int s = 0; s < nS; s++) {
                                Kval -= Dp[p + oP][s + oS] * ((*A2yp++) + (*C2yp++) + (*D2yp++)); 
                            }}
                            C_DAXPY(nocc,Kval,Cop[r + oR],1,Kmiyp[p + oP],1);
                        }}

                        // Bz
                        for (int p = 0; p < nP; p++) {
                        for (int r = 0; r < nR; r++) {
                            double Kval = 0.0;
                            for (int q = 0; q < nQ; q++) {
                            const double* A2zp = Azp + p * nQ * nR * nS + q * nR * nS + r * nS;
                            const double* C2zp = Czp + p * nQ * nR * nS + q * nR * nS + r * nS;
                            const double* D2zp = Dzp + p * nQ * nR * nS + q * nR * nS + r * nS;
                            for (int s = 0; s < nS; s++) {
                                Kval -= Dp[p + oP][s + oS] * ((*A2zp++) + (*C2zp++) + (*D2zp++)); 
                            }}
                            C_DAXPY(nocc,Kval,Cop[r + oR],1,Kmizp[p + oP],1);
                        }}
                    }

                    // C
                    if (aR == A) {

                        // Cx
                        for (int p = 0; p < nP; p++) {
                        for (int r = 0; r < nR; r++) {
                            double Kval = 0.0;
                            for (int q = 0; q < nQ; q++) {
                            const double* C2xp = Cxp + p * nQ * nR * nS + q * nR * nS + r * nS;
                            for (int s = 0; s < nS; s++) {
                                Kval += Dp[p + oP][s + oS] * (*C2xp++); 
                            }}
                            C_DAXPY(nocc,Kval,Cop[r + oR],1,Kmixp[p + oP],1);
                        }}

                        // Cy
                        for (int p = 0; p < nP; p++) {
                        for (int r = 0; r < nR; r++) {
                            double Kval = 0.0;
                            for (int q = 0; q < nQ; q++) {
                            const double* C2yp = Cyp + p * nQ * nR * nS + q * nR * nS + r * nS;
                            for (int s = 0; s < nS; s++) {
                                Kval += Dp[p + oP][s + oS] * (*C2yp++); 
                            }}
                            C_DAXPY(nocc,Kval,Cop[r + oR],1,Kmiyp[p + oP],1);
                        }}

                        // Cz
                        for (int p = 0; p < nP; p++) {
                        for (int r = 0; r < nR; r++) {
                            double Kval = 0.0;
                            for (int q = 0; q < nQ; q++) {
                            const double* C2zp = Czp + p * nQ * nR * nS + q * nR * nS + r * nS;
                            for (int s = 0; s < nS; s++) {
                                Kval += Dp[p + oP][s + oS] * (*C2zp++); 
                            }}
                            C_DAXPY(nocc,Kval,Cop[r + oR],1,Kmizp[p + oP],1);
                        }}
                    }

                    // D
                    if (aS == A) {

                        // Dx
                        for (int p = 0; p < nP; p++) {
                        for (int r = 0; r < nR; r++) {
                            double Kval = 0.0;
                            for (int q = 0; q < nQ; q++) {
                            const double* D2xp = Dxp + p * nQ * nR * nS + q * nR * nS + r * nS;
                            for (int s = 0; s < nS; s++) {
                                Kval += Dp[p + oP][s + oS] * (*D2xp++); 
                            }}
                            C_DAXPY(nocc,Kval,Cop[r + oR],1,Kmixp[p + oP],1);
                        }}

                        // Dy
                        for (int p = 0; p < nP; p++) {
                        for (int r = 0; r < nR; r++) {
                            double Kval = 0.0;
                            for (int q = 0; q < nQ; q++) {
                            const double* D2yp = Dyp + p * nQ * nR * nS + q * nR * nS + r * nS;
                            for (int s = 0; s < nS; s++) {
                                Kval += Dp[p + oP][s + oS] * (*D2yp++); 
                            }}
                            C_DAXPY(nocc,Kval,Cop[r + oR],1,Kmiyp[p + oP],1);
                        }}

                        // Dz
                        for (int p = 0; p < nP; p++) {
                        for (int r = 0; r < nR; r++) {
                            double Kval = 0.0;
                            for (int q = 0; q < nQ; q++) {
                            const double* D2zp = Dzp + p * nQ * nR * nS + q * nR * nS + r * nS;
                            for (int s = 0; s < nS; s++) {
                                Kval += Dp[p + oP][s + oS] * (*D2zp++); 
                            }}
                            C_DAXPY(nocc,Kval,Cop[r + oR],1,Kmizp[p + oP],1);
                        }}
                    }
            }}}}            

            // Jpi_x
            C_DGEMM('T','N',nmo,nocc,nso,1.0,Cp[0],nmo,Jmixp[0],nocc,0.0,Tpip[0],nocc);
            psio_->write(PSIF_HESS,"Jpi^A",(char*)Tpip[0],nmo * nocc * sizeof(double),next_Jpi,&next_Jpi);
            // Jpi_y
            C_DGEMM('T','N',nmo,nocc,nso,1.0,Cp[0],nmo,Jmiyp[0],nocc,0.0,Tpip[0],nocc);
            psio_->write(PSIF_HESS,"Jpi^A",(char*)Tpip[0],nmo * nocc * sizeof(double),next_Jpi,&next_Jpi);
            // Jpi_z
            C_DGEMM('T','N',nmo,nocc,nso,1.0,Cp[0],nmo,Jmizp[0],nocc,0.0,Tpip[0],nocc);
            psio_->write(PSIF_HESS,"Jpi^A",(char*)Tpip[0],nmo * nocc * sizeof(double),next_Jpi,&next_Jpi);
            // Kpi_x
            C_DGEMM('T','N',nmo,nocc,nso,1.0,Cp[0],nmo,Kmixp[0],nocc,0.0,Tpip[0],nocc);
            psio_->write(PSIF_HESS,"Kpi^A",(char*)Tpip[0],nmo * nocc * sizeof(double),next_Kpi,&next_Kpi);
            // Kpi_y
            C_DGEMM('T','N',nmo,nocc,nso,1.0,Cp[0],nmo,Kmiyp[0],nocc,0.0,Tpip[0],nocc);
            psio_->write(PSIF_HESS,"Kpi^A",(char*)Tpip[0],nmo * nocc * sizeof(double),next_Kpi,&next_Kpi);
            // Kpi_z
            C_DGEMM('T','N',nmo,nocc,nso,1.0,Cp[0],nmo,Kmizp[0],nocc,0.0,Tpip[0],nocc);
            psio_->write(PSIF_HESS,"Kpi^A",(char*)Tpip[0],nmo * nocc * sizeof(double),next_Kpi,&next_Kpi);
        }
    }

    boost::shared_ptr<JK> jk = JK::build_JK(basisset_, options_);
    size_t mem = 0.9 * memory_ / 8L;
    size_t per_A = 3L * nso * nso + 1L * nocc * nso;
    size_t max_A = (mem / 2L) / per_A;
    max_A = (max_A > 3 * natom ? 3 * natom : max_A);
    jk->set_memory(mem);

    // => J2pi/K2pi <= //
    {
        std::vector<boost::shared_ptr<Matrix> >& L = jk->C_left();  
        std::vector<boost::shared_ptr<Matrix> >& R = jk->C_right();  
        const std::vector<boost::shared_ptr<Matrix> >& J = jk->J();
        const std::vector<boost::shared_ptr<Matrix> >& K = jk->K();
        L.clear();
        R.clear();

        boost::shared_ptr<Matrix> Sii(new Matrix("Sii",nocc,nocc));
        double** Siip = Sii->pointer();
        boost::shared_ptr<Matrix> T(new Matrix("T",nso,nocc));
        double** Tp = T->pointer();
        boost::shared_ptr<Matrix> U(new Matrix("T",nmo,nocc));
        double** Up = U->pointer();

        for (int A = 0; A < 3 * natom; A++) {
            psio_address next_Jpi = psio_get_address(PSIO_ZERO,(A) * (size_t) nmo * nocc * sizeof(double));
            psio_->write(PSIF_HESS,"J2pi^A",(char*)Up[0],nmo*nocc*sizeof(double),next_Jpi,&next_Jpi);
        }

        for (int A = 0; A < 3 * natom; A++) {
            psio_address next_Kpi = psio_get_address(PSIO_ZERO,(A) * (size_t) nmo * nocc * sizeof(double));
            psio_->write(PSIF_HESS,"K2pi^A",(char*)Up[0],nmo*nocc*sizeof(double),next_Kpi,&next_Kpi);
        }

        for (int A = 0; A < max_A; A++) {
            L.push_back(Cocc);
            R.push_back(boost::shared_ptr<Matrix>(new Matrix("R",nso,nocc)));
        }
        jk->initialize();
        jk->print_header();

        for (int A = 0; A < 3 * natom; A+=max_A) {
            int nA = max_A;
            if (A + max_A >= 3 * natom) {
                nA = 3 * natom - A;
                L.resize(nA);
                R.resize(nA);
            }
            for (int a = 0; a < nA; a++) {
                psio_address next_Sii = psio_get_address(PSIO_ZERO,(A + a) * (size_t) nmo * nocc * sizeof(double));
                psio_->read(PSIF_HESS,"Spi^A",(char*)Siip[0],nocc*nocc*sizeof(double),next_Sii,&next_Sii);
                C_DGEMM('N','N',nso,nocc,nocc,1.0,Cop[0],nocc,Siip[0],nocc,0.0,R[a]->pointer()[0],nocc);
                L[a]->print();
                R[a]->print();
            }
            jk->compute();
            for (int a = 0; a < nA; a++) {
                C_DGEMM('N','N',nso,nocc,nso,1.0,J[a]->pointer()[0],nso,Cop[0],nocc,0.0,Tp[0],nocc);
                C_DGEMM('T','N',nmo,nocc,nso,1.0,Cp[0],nmo,Tp[0],nocc,0.0,Up[0],nocc);
                psio_address next_Jpi = psio_get_address(PSIO_ZERO,(A + a) * (size_t) nmo * nocc * sizeof(double));
                psio_->write(PSIF_HESS,"J2pi^A",(char*)Up[0],nmo*nocc*sizeof(double),next_Jpi,&next_Jpi);
                // Hermitivitize?
                C_DGEMM('N','N',nso,nocc,nso,1.0,K[a]->pointer()[0],nso,Cop[0],nocc,0.0,Tp[0],nocc);
                C_DGEMM('T','N',nmo,nocc,nso,1.0,Cp[0],nmo,Tp[0],nocc,0.0,Up[0],nocc);
                psio_address next_Kpi = psio_get_address(PSIO_ZERO,(A + a) * (size_t) nmo * nocc * sizeof(double));
                psio_->write(PSIF_HESS,"K2pi^A",(char*)Up[0],nmo*nocc*sizeof(double),next_Kpi,&next_Kpi);
            }
        } 
    }

    // => Fpi <= //
    {
        boost::shared_ptr<Matrix> Tpi(new Matrix("T",nmo,nocc));
        boost::shared_ptr<Matrix> Fpi(new Matrix("F",nmo,nocc));
        double** Tpip = Tpi->pointer();
        double** Fpip = Fpi->pointer();
        
        psio_address next_Tpi = PSIO_ZERO;    
        psio_address next_Vpi = PSIO_ZERO;    
        psio_address next_Jpi = PSIO_ZERO;    
        psio_address next_Kpi = PSIO_ZERO;    
        psio_address next_Fpi = PSIO_ZERO;    

        for (int A = 0; A < 3*natom; A++) {
            psio_->read(PSIF_HESS,"Tpi^A",(char*)Fpip[0],nmo * nocc * sizeof(double),next_Tpi,&next_Tpi);
            psio_->read(PSIF_HESS,"Vpi^A",(char*)Tpip[0],nmo * nocc * sizeof(double),next_Vpi,&next_Vpi);
            Fpi->add(Tpi);
            psio_->read(PSIF_HESS,"Jpi^A",(char*)Tpip[0],nmo * nocc * sizeof(double),next_Jpi,&next_Jpi);
            Tpi->scale(2.0);
            Fpi->add(Tpi);
            psio_->read(PSIF_HESS,"Kpi^A",(char*)Tpip[0],nmo * nocc * sizeof(double),next_Kpi,&next_Kpi);
            Tpi->scale(-1.0);
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
        psio_address next_J2pi = PSIO_ZERO;    
        psio_address next_K2pi = PSIO_ZERO;    
        psio_address next_Bai = PSIO_ZERO;    

        for (int A = 0; A < 3*natom; A++) {
            next_Fpi = psio_get_address(PSIO_ZERO,sizeof(double)*(A * (size_t) nmo * nocc + nocc * nocc));
            psio_->read(PSIF_HESS,"Fpi^A",(char*)Baip[0],nvir * nocc * sizeof(double),next_Fpi,&next_Fpi);
            next_Spi = psio_get_address(PSIO_ZERO,sizeof(double)*(A * (size_t) nmo * nocc + nocc * nocc));
            psio_->read(PSIF_HESS,"Spi^A",(char*)Taip[0],nvir * nocc * sizeof(double),next_Spi,&next_Spi);
            for (int a = 0; a < nvir; a++) {
                C_DAXPY(nocc,-evp[a],&Taip[0][a],nvir,&Baip[0][a],nvir);
            }
            next_J2pi = psio_get_address(PSIO_ZERO,sizeof(double)*(A * (size_t) nmo * nocc + nocc * nocc));
            psio_->read(PSIF_HESS,"J2pi^A",(char*)Taip[0],nvir * nocc * sizeof(double),next_J2pi,&next_J2pi);
            Tai->scale(2.0);
            Bai->add(Tai);
            next_K2pi = psio_get_address(PSIO_ZERO,sizeof(double)*(A * (size_t) nmo * nocc + nocc * nocc));
            psio_->read(PSIF_HESS,"K2pi^A",(char*)Taip[0],nvir * nocc * sizeof(double),next_K2pi,&next_K2pi);
            Tai->scale(-1.0);
            Bai->add(Tai);
            psio_->write(PSIF_HESS,"Bai^A",(char*)Baip[0],nvir * nocc * sizeof(double),next_Bai,&next_Bai);
        }
    }    

    // => CPHF (Uai) <= //
    {
        SharedWavefunction wfn(new Wavefunction(options_));
        wfn->shallow_copy(this);
                
        boost::shared_ptr<RCPHF> cphf(new RCPHF(wfn, options_));     
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
        double** Upip = Upi->pointer();    

        psio_address next_Spi = PSIO_ZERO;    
        psio_address next_Uai = PSIO_ZERO;    
        psio_address next_Upi = PSIO_ZERO;    

        for (int A = 0; A < 3*natom; A++) {
            next_Spi = psio_get_address(PSIO_ZERO,sizeof(double)*(A * (size_t) nmo * nocc + 0));
            psio_->read(PSIF_HESS,"Spi^A",(char*)Upip[0],nocc * nocc * sizeof(double),next_Spi,&next_Spi);
            C_DSCAL(nocc * (size_t) nocc, -0.5, Upip[0], 1);
            psio_->read(PSIF_HESS,"Uai^A",(char*)Upip[nocc],nvir * nocc * sizeof(double),next_Uai,&next_Uai);
            psio_->write(PSIF_HESS,"Upi^A",(char*)Upip[0],nmo * nocc * sizeof(double),next_Upi,&next_Upi);
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

        boost::shared_ptr<Matrix> Uii(new Matrix("Uii",nocc,nocc));
        double** Uiip = Uii->pointer();
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
                psio_address next_Uii = psio_get_address(PSIO_ZERO,(A + a) * (size_t) nmo * nocc * sizeof(double));
                psio_->read(PSIF_HESS,"Upi^A",(char*)Uiip[0],nocc*nocc*sizeof(double),next_Uii,&next_Uii);
                C_DGEMM('N','N',nso,nocc,nocc,1.0,Cop[0],nocc,Uiip[0],nocc,0.0,R[a]->pointer()[0],nocc);
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
        size_t max_a = memory / (2L * npi);
        max_a = (max_a > 3 * natom ? 3 * natom : max_a);
        
        boost::shared_ptr<Matrix> L(new Matrix("L",max_a * nmo, nocc));
        boost::shared_ptr<Matrix> R(new Matrix("R",max_a * nmo, nocc));
        double** Lp = L->pointer();
        double** Rp = R->pointer();

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
        
        // U^A U^B 
        for (int A = 0; A < 3 * natom; A+=max_a) {
            int nA = (A + max_a >= 3 * natom ? 3 * natom - A : max_a);
            psio_address nextA = psio_get_address(PSIO_ZERO, A * npi * sizeof(double));
            psio_->read(PSIF_HESS,"Upi^A",(char*)Lp[0],sizeof(double) * nA * npi,nextA,&nextA);
            for (int B = 0; B < 3 * natom; B+=max_a) {
                int nB = (B + max_a >= 3 * natom ? 3 * natom - B : max_a);
                psio_address nextB = psio_get_address(PSIO_ZERO, B * npi * sizeof(double));
                psio_->read(PSIF_HESS,"Upi^A",(char*)Rp[0],sizeof(double) * nB * npi,nextB,&nextB);
                for (int a = 0; a < nA; a++) {
                    for (int b = 0; b < nB; b++) {
                        Hp[A + a][B + b] -= 2.0 * C_DDOT(npi,Lp[0] + a * npi,1,Rp[0] + b * npi,1);
                    }    
                } 
            }
        }
        
        // S^A S^B 
        for (int A = 0; A < 3 * natom; A+=max_a) {
            int nA = (A + max_a >= 3 * natom ? 3 * natom - A : max_a);
            psio_address nextA = psio_get_address(PSIO_ZERO, A * npi * sizeof(double));
            psio_->read(PSIF_HESS,"Spi^A",(char*)Lp[0],sizeof(double) * nA * npi,nextA,&nextA);
            for (int B = 0; B < 3 * natom; B+=max_a) {
                int nB = (B + max_a >= 3 * natom ? 3 * natom - B : max_a);
                psio_address nextB = psio_get_address(PSIO_ZERO, B * npi * sizeof(double));
                psio_->read(PSIF_HESS,"Spi^A",(char*)Rp[0],sizeof(double) * nB * npi,nextB,&nextB);
                for (int a = 0; a < nA; a++) {
                    for (int b = 0; b < nB; b++) {
                        Hp[A + a][B + b] += 2.0 * C_DDOT(npi,Lp[0] + a * npi,1,Rp[0] + b * npi,1);
                    }    
                } 
            }
        }

        // U^A U^B \epsilon 
        for (int A = 0; A < 3 * natom; A+=max_a) {
            int nA = (A + max_a >= 3 * natom ? 3 * natom - A : max_a);
            psio_address nextA = psio_get_address(PSIO_ZERO, A * npi * sizeof(double));
            psio_->read(PSIF_HESS,"Upi^A",(char*)Lp[0],sizeof(double) * nA * npi,nextA,&nextA);
            for (int a = 0; a < nA; a++) {
                double* Tp = Lp[0]; 
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
                        Hp[A + a][B + b] += 2.0 * C_DDOT(npi,Lp[0] + a * npi,1,Rp[0] + b * npi,1);
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
                        Hp[A + a][B + b] += 2.0 * C_DDOT(npi,Lp[0] + a * npi,1,Rp[0] + b * npi,1);
                    }    
                } 
            }
        }
        
            
        // Full symmetrization
        for (int A = 0; A < 3 * natom; A++) {
            for (int B = 0; B < 3 * natom; B++) {
                Hp[A][B] = Hp[B][A] = Hp[A][B] + Hp[B][A];
            }
        }
    }

    psio_->close(PSIF_HESS,0);
    return response;
}


}} // Namespaces