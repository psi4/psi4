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

// t1-transformed 3-index fock matrix (using 3-index integrals from SCF)
void DFCoupledCluster::T1Fock(){
    long int o = ndoccact;
    long int v = nvirt;
    long int full = o+v+nfzc+nfzv;

    // Ca_L = C(1-t1^T)
    // Ca_R = C(1+t1)
    double * Catemp = (double*)malloc(nso*full*sizeof(double));
    C_DCOPY(nso*full,&Ca[0][0],1,Ca_L,1);
    C_DCOPY(nso*full,&Ca[0][0],1,Ca_R,1);
    C_DCOPY(nso*full,&Ca[0][0],1,Catemp,1);

    #pragma omp parallel for schedule (static)
    for (long int mu = 0; mu < nso; mu++) {
        for (long int a = 0; a < v; a++) {
            double dum = 0.0;
            for (long int i = 0; i < o; i++) {
                dum += Catemp[mu*full+i+nfzc] * t1[a*o+i];
            }
            Ca_L[mu*full + a + ndocc] -= dum;
        }
    }
    #pragma omp parallel for schedule (static)
    for (long int mu = 0; mu < nso; mu++) {
        for (long int i = 0; i < o; i++) {
            double dum = 0.0;
            for (long int a = 0; a < v; a++) {
                dum += Catemp[mu*full+a+ndocc] * t1[a*o+i];
            }
            Ca_R[mu*full + i + nfzc] += dum;
        }
    }
    free(Catemp);

    // (Q|rs)
    std::shared_ptr<PSIO> psio(new PSIO());
    psio->open(PSIF_DCC_QSO,PSIO_OPEN_OLD);
    psio_address addr1  = PSIO_ZERO;
    psio_address addr2  = PSIO_ZERO;
    psio_address addroo = PSIO_ZERO;
    psio_address addrov = PSIO_ZERO;
    psio_address addrvo = PSIO_ZERO;
    psio_address addrvv = PSIO_ZERO;

    long int nrows = 1;
    long int rowsize = nQ_scf;
    while ( rowsize*nso*nso > o*o*v*v ) {
        nrows++;
        rowsize = nQ_scf / nrows;
        if (nrows * rowsize < nQ_scf) rowsize++;
        if (rowsize == 1) break;
    }
    long int lastrowsize = nQ_scf - (nrows - 1L) * rowsize;
    long int * rowdims = new long int [nrows];
    for (long int i = 0; i < nrows-1; i++) rowdims[i] = rowsize;
    rowdims[nrows-1] = lastrowsize;
    for (long int row = 0; row < nrows; row++) {
        psio->read(PSIF_DCC_QSO,"Qso SCF",(char*)&integrals[0],rowdims[row]*nso*nso*sizeof(double),addr1,&addr1);
        F_DGEMM('n','n',full,nso*rowdims[row],nso,1.0,Ca_L,full,integrals,nso,0.0,tempv,full);
        for (long int q = 0; q < rowdims[row]; q++) {
            for (long int mu = 0; mu < nso; mu++) {
                C_DCOPY(full,tempv+q*nso*full+mu*full,1,integrals+q*nso*full+mu,nso);
            }
        }
        F_DGEMM('n','n',full,full*rowdims[row],nso,1.0,Ca_R,full,integrals,nso,0.0,tempv,full);
        // full Qmo
        psio->write(PSIF_DCC_QSO,"Qmo SCF",(char*)&tempv[0],rowdims[row]*full*full*sizeof(double),addr2,&addr2);
    }
    delete[] rowdims;

    // build Fock matrix
    memset((void*)Fij,'\0',o*o*sizeof(double));
    memset((void*)Fia,'\0',o*v*sizeof(double));
    memset((void*)Fai,'\0',o*v*sizeof(double));
    memset((void*)Fab,'\0',v*v*sizeof(double));

    // transform H
    double ** hp = H->pointer();
    double * h = (double*)malloc(nmo*nmo*sizeof(double));
    for (long int mu = 0; mu < nso; mu++) {
        for (long int p = 0; p < nmo; p++) {
            double dum = 0.0;
            for (long int nu = 0; nu < nso; nu++) {
                dum += Ca_L[nu*full + p + nfzc] * hp[nu][mu];
            }
            integrals[p*nso+mu] = dum;
        }
    }
    for (long int p = 0; p < nmo; p++) {
        for (long int q = 0; q < nmo; q++) {
            double dum = 0.0;
            for (long int nu = 0; nu < nso; nu++) {
                dum += Ca_R[nu*full+q+nfzc] * integrals[p*nso+nu];
            }
            h[p*nmo+q] = dum;
        }
    }

    double * temp3 = (double*)malloc(full*full*sizeof(double));

    memset((void*)temp3,'\0',full*full*sizeof(double));
    psio_address addr = PSIO_ZERO;

    nrows = 1;
    rowsize = nQ_scf;
    while ( rowsize*full*full > o*o*v*v ) {
        nrows++;
        rowsize = nQ_scf / nrows;
        if (nrows * rowsize < nQ_scf) rowsize++;
        if (rowsize == 1) break;
    }
    lastrowsize = nQ_scf - (nrows - 1L) * rowsize;
    rowdims = new long int [nrows];
    for (long int i = 0; i < nrows-1; i++) rowdims[i] = rowsize;
    rowdims[nrows-1] = lastrowsize;
    for (long int row = 0; row < nrows; row++) {
        psio->read(PSIF_DCC_QSO,"Qmo SCF",(char*)&integrals[0],rowdims[row]*full*full*sizeof(double),addr,&addr);
        for (long int q = 0; q < rowdims[row]; q++) {
            // sum k (q|rk) (q|ks)
            F_DGEMM('n','n',full,full,ndocc,-1.0,integrals+q*full*full,full,integrals+q*full*full,full,1.0,temp3,full);

            // sum k (q|kk) (q|rs)
            double dum = 0.0;
            for (long int k = 0; k < ndocc; k++) {
                dum += integrals[q*full*full+k*full + k];
            }
            C_DAXPY(full*full,2.0 * dum,integrals+q*full*full,1,temp3,1);
        }
    }
    delete[] rowdims;
    psio->close(PSIF_DCC_QSO,1);

    // Fij
    for (long int i = 0; i < o; i++) {
        for (long int j = 0; j < o; j++) {
            Fij[i*o+j] = h[i*nmo+j] + temp3[(i+nfzc)*full+(j+nfzc)];
        }
    }

    // Fia
    for (long int i = 0; i < o; i++) {
        for (long int a = 0; a < v; a++) {
            Fia[i*v+a] = h[i*nmo+a+o] + temp3[(i+nfzc)*full+(a+ndocc)];
        }
    }

    // Fai
    for (long int a = 0; a < v; a++) {
        for (long int i = 0; i < o; i++) {
            Fai[a*o+i] = h[(a+o)*nmo+i] + temp3[(a+ndocc)*full+(i+nfzc)];
        }
    }

    // Fab
    for (long int a = 0; a < v; a++) {
        for (long int b = 0; b < v; b++) {
            Fab[a*v+b] = h[(a+o)*nmo+b+o] + temp3[(a+ndocc)*full+(b+ndocc)];
        }
    }

    // replace eps
    for (long int i = 0; i < o; i++) {
        eps[i] = Fij[i*o+i];
    }
    for (long int a = 0; a < v; a++) {
        eps[a+o] = Fab[a*v+a];
    }

    free(h);
    free(temp3);
}

void DFCoupledCluster::T1Integrals(){
    long int o = ndoccact;
    long int v = nvirt;
    long int full = o+v+nfzc+nfzv;

    // Ca_L = C(1-t1^T)
    // Ca_R = C(1+t1)
    double * Catemp = (double*)malloc(nso*full*sizeof(double));
    C_DCOPY(nso*full,&Ca[0][0],1,Ca_L,1);
    C_DCOPY(nso*full,&Ca[0][0],1,Ca_R,1);
    C_DCOPY(nso*full,&Ca[0][0],1,Catemp,1);

    #pragma omp parallel for schedule (static)
    for (long int mu = 0; mu < nso; mu++) {
        for (long int a = 0; a < v; a++) {
            double dum = 0.0;
            for (long int i = 0; i < o; i++) {
                dum += Catemp[mu*full+i+nfzc] * t1[a*o+i];
            }
            Ca_L[mu*full + a + ndocc] -= dum;
        }
    }
    #pragma omp parallel for schedule (static)
    for (long int mu = 0; mu < nso; mu++) {
        for (long int i = 0; i < o; i++) {
            double dum = 0.0;
            for (long int a = 0; a < v; a++) {
                dum += Catemp[mu*full+a+ndocc] * t1[a*o+i];
            }
            Ca_R[mu*full + i + nfzc] += dum;
        }
    }
    free(Catemp);

    // (Q|rs)
    std::shared_ptr<PSIO> psio(new PSIO());
    psio->open(PSIF_DCC_QSO,PSIO_OPEN_OLD);
    psio_address addr1  = PSIO_ZERO;
    psio_address addrvo = PSIO_ZERO;
    long int nrows = 1;
    long int rowsize = nQ;
    while ( rowsize*nso*nso > o*o*v*v ) {
        nrows++;
        rowsize = nQ / nrows;
        if (nrows * rowsize < nQ) rowsize++;
        if (rowsize == 1) break;
    }
    long int lastrowsize = nQ - (nrows - 1L) * rowsize;
    long int * rowdims = new long int [nrows];
    for (long int i = 0; i < nrows-1; i++) rowdims[i] = rowsize;
    rowdims[nrows-1] = lastrowsize;
    for (long int row = 0; row < nrows; row++) {
        psio->read(PSIF_DCC_QSO,"Qso CC",(char*)&integrals[0],rowdims[row]*nso*nso*sizeof(double),addr1,&addr1);
        F_DGEMM('n','n',full,nso*rowdims[row],nso,1.0,Ca_L,full,integrals,nso,0.0,tempv,full);
        for (long int q = 0; q < rowdims[row]; q++) {
            for (long int mu = 0; mu < nso; mu++) {
                C_DCOPY(full,tempv+q*nso*full+mu*full,1,integrals+q*nso*full+mu,nso);
            }
        }
        F_DGEMM('n','n',full,full*rowdims[row],nso,1.0,Ca_R,full,integrals,nso,0.0,tempv,full);

        // Qoo
        #pragma omp parallel for schedule (static)
        for (long int q = 0; q < rowdims[row]; q++) {
            for (long int i = 0; i < o; i++) {
                for (long int j = 0; j < o; j++) {
                    Qoo[(q+rowdims[0]*row)*o*o+i*o+j] = tempv[q*full*full+(i+nfzc)*full+(j+nfzc)];
                }
            }
        }
        // Qov
        #pragma omp parallel for schedule (static)
        for (long int q = 0; q < rowdims[row]; q++) {
            for (long int i = 0; i < o; i++) {
                for (long int a = 0; a < v; a++) {
                    Qov[(q+rowdims[0]*row)*o*v+i*v+a] = tempv[q*full*full+(i+nfzc)*full+(a+ndocc)];
                }
            }
        }
        // Qvo
        #pragma omp parallel for schedule (static)
        for (long int q = 0; q < rowdims[row]; q++) {
            for (long int a = 0; a < v; a++) {
                for (long int i = 0; i < o; i++) {
                    integrals[q*o*v+a*o+i] = tempv[q*full*full+(a+ndocc)*full+(i+nfzc)];
                }
            }
        }
        psio->write(PSIF_DCC_QSO,"qvo",(char*)&integrals[0],rowdims[row]*o*v*sizeof(double),addrvo,&addrvo);
        // Qvv
        #pragma omp parallel for schedule (static)
        for (long int q = 0; q < rowdims[row]; q++) {
            for (long int a = 0; a < v; a++) {
                for (long int b = 0; b < v; b++) {
                    Qvv[(q+rowdims[0]*row)*v*v+a*v+b] = tempv[q*full*full+(a+ndocc)*full+(b+ndocc)];
                }
            }
        }
    }
    delete[] rowdims;
    psio->close(PSIF_DCC_QSO,1);
}

}}
