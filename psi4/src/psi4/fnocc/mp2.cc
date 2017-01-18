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
#include "psi4/libtrans/mospace.h"
#include "psi4/libtrans/integraltransform.h"
#include "psi4/libiwl/iwl.h"
#include "psi4/psifiles.h"
#ifdef _OPENMP
    #include<omp.h>
#else
    #define omp_get_wtime() 0.0
#endif

#include"blas.h"
#include"ccsd.h"

using namespace psi;
namespace psi{ namespace fnocc{
void SortOVOV(struct iwlbuf *Buf,int nfzc,int nfzv,int norbs,int ndoccact,int nvirt);
}}

namespace psi{ namespace fnocc{

/**
  * canonical MP2
  */
void CoupledCluster::MP2(){
    int o = ndoccact;
    int v = nvirt;

    // transform integrals
    outfile->Printf("\n");
    outfile->Printf("        ==> Transform (OV|OV) integrals <==\n");
    outfile->Printf("\n");
    std::shared_ptr<psi::Wavefunction> wfn = reference_wavefunction_;
    std::vector<std::shared_ptr<MOSpace> > spaces;
    spaces.push_back(MOSpace::occ);
    spaces.push_back(MOSpace::vir);
    std::shared_ptr<IntegralTransform>
        ints(new IntegralTransform(wfn,
                                   spaces,
                                   IntegralTransform::Restricted,
                                   IntegralTransform::IWLOnly,
                                   IntegralTransform::QTOrder,
                                   IntegralTransform::None,
                                   false));
    ints->set_keep_dpd_so_ints(1);
    ints->set_keep_iwl_so_ints(1);
    ints->initialize();
    ints->transform_tei(MOSpace::occ, MOSpace::vir, MOSpace::occ, MOSpace::vir);
    ints.reset();

    // sort integrals
    outfile->Printf("\n");
    outfile->Printf("        ==> Sort (OV|OV) integrals <==\n");
    outfile->Printf("\n");
    struct iwlbuf Buf;
    iwl_buf_init(&Buf,PSIF_MO_TEI,0.0,1,1);
    SortOVOV(&Buf,nfzc,nfzv,nfzc+nfzv+ndoccact+nvirt,ndoccact,nvirt);
    iwl_buf_close(&Buf,1);


    // energy
    double * v2 = (double*)malloc(o*o*v*v*sizeof(double));

    std::shared_ptr<PSIO> psio(new PSIO());
    psio->open(PSIF_DCC_IAJB,PSIO_OPEN_OLD);
    psio->read_entry(PSIF_DCC_IAJB,"E2iajb",(char*)&v2[0],o*o*v*v*sizeof(double));
    psio->close(PSIF_DCC_IAJB,0);

    emp2_os = emp2_ss = 0.0;
    for (int i = 0; i < o; i++) {
        double di = eps[i];
        for (int j = 0; j < o; j++) {
            double dij = di + eps[j];
            for (int a = 0; a < v; a++) {
                double dija = dij - eps[a+o];
                for (int b = 0; b < v; b++) {
                    double dijab = dija - eps[b+o];
                    double val1 = v2[i*o*v*v+a*o*v+j*v+b];
                    emp2_os += val1 * val1 / dijab;
                    emp2_ss += (val1 - v2[i*o*v*v+b*o*v+j*v+a]) * val1 / dijab;
                }
            }
        }
    }
    emp2 = emp2_os + emp2_ss;

    outfile->Printf("        OS MP2 correlation energy:       %20.12lf\n",emp2_os);
    outfile->Printf("        SS MP2 correlation energy:       %20.12lf\n",emp2_ss);
    outfile->Printf("        MP2 correlation energy:          %20.12lf\n",emp2);
    outfile->Printf("      * MP2 total energy:                %20.12lf\n",emp2+escf);
    outfile->Printf("\n");

    free(v2);
}

}} // end of namespaces
