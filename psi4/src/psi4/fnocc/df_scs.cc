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

void DFCoupledCluster::SCS_CCSD(){

    long int v  = nvirt;
    long int o  = ndoccact;
    long int rs = nmo;

    double ssenergy = 0.0;
    double osenergy = 0.0;

    // df (ia|bj) formerly E2klcd
    F_DGEMM('n','t',o*v,o*v,nQ,1.0,Qov,o*v,Qov,o*v,0.0,integrals,o*v);

    if (t2_on_disk){
        std::shared_ptr<PSIO> psio (new PSIO());
        psio->open(PSIF_DCC_T2,PSIO_OPEN_OLD);
        psio->read_entry(PSIF_DCC_T2,"t2",(char*)&tempv[0],o*o*v*v*sizeof(double));
        psio->close(PSIF_DCC_T2,1);
        tb = tempv;
    }

    for (long int a = o; a < rs; a++){
        for (long int b = o; b < rs; b++){
            for (long int i = 0; i < o; i++){
                for (long int j = 0; j < o; j++){

                    long int ijab = (a-o)*v*o*o+(b-o)*o*o+i*o+j;
                    long int iajb = i*v*v*o+(a-o)*v*o+j*v+(b-o);
                    long int jaib = iajb + (i-j)*v*(1-v*o);

                    osenergy += integrals[iajb]*(tb[ijab]+t1[(a-o)*o+i]*t1[(b-o)*o+j]);
                    ssenergy += integrals[iajb]*(tb[ijab]-tb[(b-o)*o*o*v+(a-o)*o*o+i*o+j]);
                    ssenergy += integrals[iajb]*(t1[(a-o)*o+i]*t1[(b-o)*o+j]-t1[(b-o)*o+i]*t1[(a-o)*o+j]);
                }
            }
        }
    }
    eccsd_os = osenergy;
    eccsd_ss = ssenergy;
    eccsd    = eccsd_os + eccsd_ss;
}

void DFCoupledCluster::SCS_MP2(){

    long int v  = nvirt;
    long int o  = ndoccact;
    long int rs = nmo;

    double ssenergy = 0.0;
    double osenergy = 0.0;

    // df (ia|bj) formerly E2klcd
    F_DGEMM('n','t',o*v,o*v,nQ,1.0,Qov,o*v,Qov,o*v,0.0,integrals,o*v);

    if (t2_on_disk){
        std::shared_ptr<PSIO> psio (new PSIO());
        psio->open(PSIF_DCC_T2,PSIO_OPEN_OLD);
        psio->read_entry(PSIF_DCC_T2,"t2",(char*)&tempv[0],o*o*v*v*sizeof(double));
        psio->close(PSIF_DCC_T2,1);
        tb = tempv;
    }

    for (long int a = o; a < rs; a++){
        for (long int b = o; b < rs; b++){
            for (long int i = 0; i < o; i++){
                for (long int j = 0; j < o; j++){

                    long int ijab = (a-o)*v*o*o+(b-o)*o*o+i*o+j;
                    long int iajb = i*v*v*o+(a-o)*v*o+j*v+(b-o);
                    long int jaib = iajb + (i-j)*v*(1-v*o);
                    osenergy += integrals[iajb]*tb[ijab];
                    ssenergy += integrals[iajb]*(tb[ijab]-tb[(b-o)*o*o*v+(a-o)*o*o+i*o+j]);
                }
            }
        }
    }
    emp2_os = osenergy;
    emp2_ss = ssenergy;
    emp2    = emp2_os + emp2_ss;
}

}}
