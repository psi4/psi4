/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
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

#include "psi4/psi4-dec.h"
#include "psi4/libmints/vector.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/wavefunction.h"
#include"psi4/libqt/qt.h"
#include<sys/times.h>
#include "psi4/libciomr/libciomr.h"
#ifdef _OPENMP
    #include<omp.h>
#else
    #define omp_get_wtime() 0.0
#endif

#include"blas.h"
#include"ccsd.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/basisset_parser.h"
#include "psi4/lib3index/3index.h"

using namespace psi;


namespace psi{ namespace fnocc{

void DFCoupledCluster::CCResidual(){
    bool timer = options_.get_bool("CC_TIMINGS");
    long int o = ndoccact;
    long int v = nvirt;

    double start;

    std::shared_ptr<PSIO> psio (new PSIO());

    // C2 = -1/2 t(bc,kj) [ (ki|ac) - 1/2 t(ad,li) (kd|lc) ]
    //      +    t(bc,ki) [ (kj|ac) - 1/2 t(ad,lj) (kd|lc) ]
    if (timer) start = omp_get_wtime();
    F_DGEMM('n','t',o*v,o*v,nQ,1.0,Qov,o*v,Qov,o*v,0.0,integrals,o*v);
    if (t2_on_disk){
        psio->open(PSIF_DCC_T2,PSIO_OPEN_OLD);
        psio->read_entry(PSIF_DCC_T2,"t2",(char*)&tempv[0],o*o*v*v*sizeof(double));
        psio->close(PSIF_DCC_T2,1);
        tb = tempv;
    }
    #pragma omp parallel for schedule (static)
    for (int a = 0; a < v; a++) {
        for (int i = 0; i < o; i++) {
            for (int l = 0; l < o; l++) {
                for (int d = 0; d < v; d++) {
                    tempt[a*o*o*v+i*o*v+l*v+d] = tb[a*o*o*v+d*o*o+l*o+i];
                }
            }
        }
    }
    #pragma omp parallel for schedule (static)
    for (int l = 0; l < o; l++) {
        for (int d = 0; d < v; d++) {
            for (int k = 0; k < o; k++) {
                for (int c = 0; c < v; c++) {
                    tempv[l*o*v*v+d*o*v+k*v+c] = integrals[k*o*v*v+d*o*v+l*v+c];
                }
            }
        }
    }
    F_DGEMM('n','n',o*v,o*v,o*v,-0.5,tempv,o*v,tempt,o*v,0.0,integrals,o*v);
    F_DGEMM('n','t',v*v,o*o,nQ,1.0,Qvv,v*v,Qoo,o*o,0.0,tempv,v*v);
    #pragma omp parallel for schedule (static)
    for (int a = 0; a < v; a++) {
        for (int i = 0; i < o; i++) {
            for (int k = 0; k < o; k++) {
                for (int c = 0; c < v; c++) {
                    integrals[a*o*o*v+i*o*v+k*v+c] += tempv[k*o*v*v+i*v*v+a*v+c];
                }
            }
        }
    }
    if (t2_on_disk){
        psio->open(PSIF_DCC_T2,PSIO_OPEN_OLD);
        psio->read_entry(PSIF_DCC_T2,"t2",(char*)&tempv[0],o*o*v*v*sizeof(double));
        psio->close(PSIF_DCC_T2,1);
        tb = tempv;
    }
    #pragma omp parallel for schedule (static)
    for (int b = 0; b < v; b++) {
        for (int j = 0; j < o; j++) {
            for (int k = 0; k < o; k++) {
                for (int c = 0; c < v; c++) {
                    tempt[b*o*o*v+j*o*v+k*v+c] = tb[b*o*o*v+c*o*o+k*o+j];
                }
            }
        }
    }
    F_DGEMM('t','n',o*v,o*v,o*v,-1.0,integrals,o*v,tempt,o*v,0.0,tempv,o*v);
    #pragma omp parallel for schedule (static)
    for (int a = 0; a < v; a++) {
        for (int b = 0; b < v; b++) {
            for (int i = 0; i < o; i++) {
                for (int j = 0; j < o; j++) {
                    tempt[a*o*o*v+b*o*o+i*o+j] = 0.5 * tempv[b*o*o*v+j*o*v+a*o+i] + tempv[b*o*o*v+i*o*v+a*o+j];
                }
            }
        }
    }

    // first contribution to residual
    psio->open(PSIF_DCC_R2,PSIO_OPEN_NEW);
    psio->write_entry(PSIF_DCC_R2,"residual",(char*)&tempt[0],o*o*v*v*sizeof(double));
    psio->close(PSIF_DCC_R2,1);
    if (timer) {
        outfile->Printf("\n");
        outfile->Printf("        C2 = -1/2 t(b,c,k,j) [ (ki|ac) - 1/2 t(a,d,l,i) (kd|lc) ]\n");
        outfile->Printf("                + t(b,c,k,i) [ (kj|ac) - 1/2 t(a,d,l,j) (kd|lc) ]       %6.2lf\n",omp_get_wtime()-start);
        start = omp_get_wtime();
    }

    // D2: 1/2 U(b,c,j,k) [ L(a,i,k,c) + 1/2 U(a,d,i,l) L(l,d,k,c) ]
    F_DGEMM('n','t',o*v,o*v,nQ,1.0,Qov,o*v,Qov,o*v,0.0,integrals,o*v);
    C_DCOPY(o*o*v*v,integrals,1,tempv,1);
    #pragma omp parallel for schedule (static)
    for (int l = 0; l < o; l++) {
        for (int d = 0; d < v; d++) {
            for (int k = 0; k < o; k++) {
                for (int c = 0; c < v; c++) {
                    tempv[l*o*v*v+d*o*v+k*v+c] -= 0.5 * integrals[l*o*v*v+c*o*v+k*v+d];
                }
            }
        }
    }
    if (t2_on_disk){
        psio->open(PSIF_DCC_T2,PSIO_OPEN_OLD);
        psio->read_entry(PSIF_DCC_T2,"t2",(char*)&integrals[0],o*o*v*v*sizeof(double));
        psio->close(PSIF_DCC_T2,1);
        tb = integrals;
    }
    #pragma omp parallel for schedule (static)
    for (int l = 0; l < o; l++) {
        for (int d = 0; d < v; d++) {
            for (int a = 0; a < v; a++) {
                for (int i = 0; i < o; i++) {
                    tempt[l*o*v*v+d*o*v+a*o+i] = 2.0 * tb[a*o*o*v+d*o*o+i*o+l]-tb[a*o*o*v+d*o*o+l*o+i];
                }
            }
        }
    }
    F_DGEMM('n','t',o*v,o*v,o*v,1.0,tempv,o*v,tempt,o*v,0.0,integrals,o*v);
    psio->open(PSIF_DCC_QSO,PSIO_OPEN_OLD);
    psio->read_entry(PSIF_DCC_QSO,"qvo",(char*)&tempv[0],nQ*o*v*sizeof(double));
    psio->close(PSIF_DCC_QSO,1);
    F_DGEMM('n','t',o*v,o*v,nQ,2.0,Qov,o*v,tempv,o*v,1.0,integrals,o*v);
    F_DGEMM('n','t',o*o,v*v,nQ,-1.0,Qoo,o*o,Qvv,v*v,0.0,tempv,o*o);
    #pragma omp parallel for schedule (static)
    for (int a = 0; a < v; a++) {
        for (int i = 0; i < o; i++) {
            for (int k = 0; k < o; k++) {
                for (int c = 0; c < v; c++) {
                    integrals[a*o*o*v+i*o*v+k*v+c] += tempv[a*o*o*v+c*o*o+k*o+i];
                }
            }
        }
    }
    if (t2_on_disk){
        psio->open(PSIF_DCC_T2,PSIO_OPEN_OLD);
        psio->read_entry(PSIF_DCC_T2,"t2",(char*)&tempv[0],o*o*v*v*sizeof(double));
        psio->close(PSIF_DCC_T2,1);
        tb = tempv;
    }
    #pragma omp parallel for schedule (static)
    for (int k = 0; k < o; k++) {
        for (int c = 0; c < v; c++) {
            for (int b = 0; b < v; b++) {
                for (int j = 0; j < o; j++) {
                    tempt[k*o*v*v+c*o*v+b*o+j] = 2.0 * tb[b*o*o*v+c*o*o+j*o+k] - tb[b*o*o*v+c*o*o+k*o+j];
                }
            }
        }
    }
    F_DGEMM('n','n',o*v,o*v,o*v,0.5,tempt,o*v,integrals,o*v,0.0,tempv,o*v);
    #pragma omp parallel for schedule (static)
    for (int a = 0; a < v; a++) {
        for (int b = 0; b < v; b++) {
            for (int i = 0; i < o; i++) {
                for (int j = 0; j < o; j++) {
                    tempt[a*o*o*v+b*o*o+i*o+j] = tempv[a*o*o*v+i*o*v+b*o+j];
                }
            }
        }
    }
    psio->open(PSIF_DCC_R2,PSIO_OPEN_OLD);
    psio->read_entry(PSIF_DCC_R2,"residual",(char*)&tempv[0],o*o*v*v*sizeof(double));
    C_DAXPY(o*o*v*v,1.0,tempv,1,tempt,1);
    psio->write_entry(PSIF_DCC_R2,"residual",(char*)&tempt[0],o*o*v*v*sizeof(double));
    psio->close(PSIF_DCC_R2,1);
    if (timer) {
        outfile->Printf("        D2 =  1/2 U(b,c,j,k) [ L(a,i,k,c) + 1/2 U(a,d,i,l) L(l,d,k,c) ] %6.2lf\n",omp_get_wtime()-start);
        start = omp_get_wtime();
    }

    // E2 a: t(ac,ij) [ F(bc) - U(bd,kl) (ld|kc) ]
    if (t2_on_disk){
        psio->open(PSIF_DCC_T2,PSIO_OPEN_OLD);
        psio->read_entry(PSIF_DCC_T2,"t2",(char*)&tempv[0],o*o*v*v*sizeof(double));
        psio->close(PSIF_DCC_T2,1);
        tb = tempv;
    }
    C_DCOPY(o*o*v*v,tb,1,tempt,1);
    #pragma omp parallel for schedule (static)
    for (int b = 0; b < v; b++) {
        for (int d = 0; d < v; d++) {
            for (int k = 0; k < o; k++) {
                C_DAXPY(o,-0.5,tb+b*o*o*v+d*o*o+k,o,tempt+b*o*o*v+d*o*o+k*o,1);
            }
        }
    }
    F_DGEMM('n','t',o*v,o*v,nQ,1.0,Qov,o*v,Qov,o*v,0.0,integrals,o*v);
    #pragma omp parallel for schedule (static)
    for (int c = 0; c < v; c++) {
        for (int d = 0; d < v; d++) {
            for (int k = 0; k < o; k++) {
                for (int l = 0; l < o; l++) {
                    tempv[c*o*o*v+d*o*o+k*o+l] = integrals[l*o*v*v+d*o*v+k*v+c];
                }
            }
        }
    }
    // overwriting Fab here, but it gets rebuilt every iteration anyway.
    F_DGEMM('t','n',v,v,o*o*v,-2.0,tempv,o*o*v,tempt,o*o*v,1.0,Fab,v);
    if (t2_on_disk){
        psio->open(PSIF_DCC_T2,PSIO_OPEN_OLD);
        psio->read_entry(PSIF_DCC_T2,"t2",(char*)&tempv[0],o*o*v*v*sizeof(double));
        psio->close(PSIF_DCC_T2,1);
        tb = tempv;
    }
    #pragma omp parallel for schedule (static)
    for (int c = 0; c < v; c++) {
        for (int a = 0; a < v; a++) {
            for (int i = 0; i < o; i++) {
                for (int j = 0; j < o; j++) {
                    tempt[c*o*o*v+a*o*o+i*o+j] = tb[a*o*o*v+c*o*o+i*o+j];
                }
            }
        }
    }
    F_DGEMM('n','n',o*o*v,v,v,1.0,tempt,o*o*v,Fab,v,0.0,tempv,o*o*v);
    #pragma omp parallel for schedule (static)
    for (int a = 0; a < v; a++) {
        for (int b = 0; b < v; b++) {
            for (int i = 0; i < o; i++) {
                for (int j = 0; j < o; j++) {
                    tempt[a*o*o*v+b*o*o+i*o+j] = tempv[b*o*o*v+a*o*o+i*o+j];
                }
            }
        }
    }
    psio->open(PSIF_DCC_R2,PSIO_OPEN_OLD);
    psio->read_entry(PSIF_DCC_R2,"residual",(char*)&tempv[0],o*o*v*v*sizeof(double));
    C_DAXPY(o*o*v*v,1.0,tempv,1,tempt,1);
    psio->write_entry(PSIF_DCC_R2,"residual",(char*)&tempt[0],o*o*v*v*sizeof(double));
    psio->close(PSIF_DCC_R2,1);
    if (timer) {
        outfile->Printf("        E2 =      t(a,c,i,j) [ F(b,c) - U(b,d,k,l) (ld|kc) ]            %6.2lf\n",omp_get_wtime()-start);
        start = omp_get_wtime();
    }

    // E2 b: -t(a,b,i,k) [ F(kj) - U(c,d,l,j) (kd|lc) ]
    // note that (kd|lc) should still be in integrals buffer
    if (t2_on_disk){
        psio->open(PSIF_DCC_T2,PSIO_OPEN_OLD);
        psio->read_entry(PSIF_DCC_T2,"t2",(char*)&tempv[0],o*o*v*v*sizeof(double));
        psio->close(PSIF_DCC_T2,1);
        tb = tempv;
    }
    #pragma omp parallel for schedule (static)
    for (int j = 0; j < o; j++) {
        for (int d = 0; d < v; d++) {
            for (int l = 0; l < o; l++) {
                for (int c = 0; c < v; c++) {
                    tempt[j*o*v*v+d*o*v+l*v+c] = (2.0 * tb[c*o*o*v+d*o*o+l*o+j] - tb[c*o*o*v+d*o*o+j*o+l] );
                }
            }
        }
    }
    // overwriting Fij here, but it gets rebuilt every iteration anyway.
    F_DGEMM('t','n',o,o,o*v*v,1.0,tempt,o*v*v,integrals,o*v*v,1.0,Fij,o);

    psio->open(PSIF_DCC_R2,PSIO_OPEN_OLD);
    psio->read_entry(PSIF_DCC_R2,"residual",(char*)&tempt[0],o*o*v*v*sizeof(double));
    F_DGEMM('n','n',o,o*v*v,o,-1.0,Fij,o,tb,o,1.0,tempt,o);

    // R2 = R2 + P(ia,jb) R2
    C_DCOPY(o*o*v*v,tempt,1,integrals,1);
    #pragma omp parallel for schedule (static)
    for (int a = 0; a < v; a++) {
        for (int b = 0; b < v; b++) {
            for (int i = 0; i < o; i++) {
                for (int j = 0; j < o; j++) {
                    integrals[a*o*o*v+b*o*o+i*o+j] += tempt[b*o*o*v+a*o*o+j*o+i];
                }
            }
        }
    }
    psio->write_entry(PSIF_DCC_R2,"residual",(char*)&integrals[0],o*o*v*v*sizeof(double));
    psio->close(PSIF_DCC_R2,1);
    if (timer) {
        outfile->Printf("                - t(a,b,i,k) [ F(k,j) - U(c,d,l,j) (kd|lc) ]            %6.2lf\n",omp_get_wtime()-start);
        start = omp_get_wtime();
    }

    // B2 = t(ab,kl) [ (ki|lj) + t(cd,ij) (kc|ld) ]
    F_DGEMM('n','t',o*v,o*v,nQ,1.0,Qov,o*v,Qov,o*v,0.0,integrals,o*v);
    #pragma omp parallel for schedule (static)
    for (int k = 0; k < o; k++) {
        for (int l = 0; l < o; l++) {
            for (int c = 0; c < v; c++) {
                for (int d = 0; d < v; d++) {
                    tempv[k*o*v*v+l*v*v+c*v+d] = integrals[k*v*v*o+c*o*v+l*v+d];
                }
            }
        }
    }
    F_DGEMM('n','t',o*o,o*o,nQ,1.0,Qoo,o*o,Qoo,o*o,0.0,integrals,o*o);
    #pragma omp parallel for schedule (static)
    for (int k = 0; k < o; k++) {
        for (int i = 0; i < o; i++) {
            for (int l = 0; l < o; l++) {
                for (int j = 0; j < o; j++) {
                    tempt[k*o*o*o+l*o*o+i*o+j] = integrals[k*o*o*o+i*o*o+l*o+j];
                }
            }
        }
    }
    if (t2_on_disk){
        psio->open(PSIF_DCC_T2,PSIO_OPEN_OLD);
        psio->read_entry(PSIF_DCC_T2,"t2",(char*)&integrals[0],o*o*v*v*sizeof(double));
        psio->close(PSIF_DCC_T2,1);
        tb = integrals;
    }
    F_DGEMM('n','n',o*o,o*o,v*v,1.0,tb,o*o,tempv,v*v,1.0,tempt,o*o);
    if (t2_on_disk){
        psio->open(PSIF_DCC_T2,PSIO_OPEN_OLD);
        psio->read_entry(PSIF_DCC_T2,"t2",(char*)&tempv[0],o*o*v*v*sizeof(double));
        psio->close(PSIF_DCC_T2,1);
        tb = tempv;
    }
    F_DGEMM('n','n',o*o,v*v,o*o,1.0,tempt,o*o,tb,o*o,0.0,integrals,o*o);

    psio->open(PSIF_DCC_R2,PSIO_OPEN_OLD);
    psio->read_entry(PSIF_DCC_R2,"residual",(char*)&tempt[0],o*o*v*v*sizeof(double));
    C_DAXPY(o*o*v*v,1.0,tempt,1,integrals,1);
    psio->write_entry(PSIF_DCC_R2,"residual",(char*)&integrals[0],o*o*v*v*sizeof(double));
    psio->close(PSIF_DCC_R2,1);

    if (timer) {
        outfile->Printf("        B2 =      t(a,b,k,l) [ (ki|lj) + t(c,d,i,j) (kc|ld) ]           %6.2lf\n",omp_get_wtime()-start);
        start = omp_get_wtime();
    }

    // now singles residual:

    // D1: F(ai)
    C_DCOPY(o*v,Fai,1,w1,1);

    // A1 (G):  U(c,d,k,l) (ad|kc)
    #pragma omp parallel for schedule (static)
    for (int d = 0; d < v; d++) {
        for (int i = 0; i < o; i++) {
            for (int k = 0; k < o; k++) {
                for (int c = 0; c < v; c++) {
                    tempt[d*o*o*v+i*o*v+k*v+c] = (2.0*tb[c*o*o*v+d*o*o+k*o+i] - tb[c*o*o*v+d*o*o+i*o+k]);
                }
            }
        }
    }
    F_DGEMM('t','n',o*v,nQ,o*v,1.0,tempt,o*v,Qov,o*v,0.0,tempv,o*v);
    #pragma omp parallel for schedule (static)
    for (int q = 0; q < nQ; q++) {
        for (int a = 0; a < v; a++) {
            C_DCOPY(v,Qvv+q*v*v+a*v,1,integrals+q*v*v+a,v);
        }
    }
    F_DGEMM('n','t',o,v,v*nQ,1.0,tempv,o,integrals,v,1.0,w1,o);

    if (timer) {
        outfile->Printf("        A1 =      U(c,d,k,l) (ad|kc)                                    %6.2lf\n",omp_get_wtime()-start);
        start = omp_get_wtime();
    }

    // B1 (H): -U(a,c,k,l) (ki|lc)
    F_DGEMM('n','t',o*v,o*o,nQ,1.0,Qov,o*v,Qoo,o*o,0.0,integrals,o*v);
    #pragma omp parallel for schedule (static)
    for (int i = 0; i < o; i++) {
        for (int c = 0; c < v; c++) {
            for (int k = 0; k < o; k++) {
                for (int l = 0; l < o; l++) {
                    tempv[i*o*o*v+c*o*o+k*o+l] = integrals[k*o*o*v+i*o*v+l*v+c];
                }
            }
        }
    }
    if (t2_on_disk){
        psio->open(PSIF_DCC_T2,PSIO_OPEN_OLD);
        psio->read_entry(PSIF_DCC_T2,"t2",(char*)&integrals[0],o*o*v*v*sizeof(double));
        psio->close(PSIF_DCC_T2,1);
        tb = integrals;
    }
    C_DCOPY(o*o*v*v,tb,1,tempt,1);
    #pragma omp parallel for schedule (static)
    for (int a = 0; a < v; a++) {
        for (int c = 0; c < v; c++) {
            for (int k = 0; k < o; k++) {
                C_DAXPY(o,-0.5,tb+a*o*o*v+c*o*o+k,o,tempt+a*o*o*v+c*o*o+k*o,1);
            }
        }
    }
    F_DGEMM('t','n',o,v,o*o*v,-2.0,tempv,o*o*v,tempt,o*o*v,1.0,w1,o);

    if (timer) {
        outfile->Printf("        B1 =    - U(a,c,k,l) (ki|lc)                                    %6.2lf\n",omp_get_wtime()-start);
        start = omp_get_wtime();
    }

    // C1
    if (t2_on_disk){
        psio->open(PSIF_DCC_T2,PSIO_OPEN_OLD);
        psio->read_entry(PSIF_DCC_T2,"t2",(char*)&tempv[0],o*o*v*v*sizeof(double));
        psio->close(PSIF_DCC_T2,1);
        tb = tempv;
    }
    #pragma omp parallel for schedule (static)
    for (int a = 0; a < v; a++) {
        for (int i = 0; i < o; i++) {
            double dum = 0.0;
            for (int k = 0; k < o; k++) {
                for (int c = 0; c < v; c++) {
                    dum += Fia[k*v+c] * (2.0*tb[a*o*o*v+c*o*o+i*o+k] - tb[a*o*o*v+c*o*o+k*o+i]);
                }
            }
            w1[a*o+i] += dum;
        }
    }

    if (timer) {
        outfile->Printf("        C1 =      F(k,c) U(a,c,i,k)                                     %6.2lf\n",omp_get_wtime()-start);
        start = omp_get_wtime();
    }

    Vabcd1();
    if (timer) {
        outfile->Printf("        A2 =      t(c,d,i,j) (ac|bd)                                    %6.2lf\n",omp_get_wtime()-start);
    }
}

}}
