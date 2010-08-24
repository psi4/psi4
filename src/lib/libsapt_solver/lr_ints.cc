/* 
 *  LR_INTS.CC 
 *
 */

#ifdef HAVE_MKL
#include <mkl.h>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <vector>
#include <utility>

#include <psifiles.h>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>
#include <libchkpt/chkpt.hpp>
#include <libipv1/ip_lib.h>
#include <libiwl/iwl.hpp>
#include <libqt/qt.h>
#include <psifiles.h>

#include "sapt.h"
#include "structs.h"

#include <libmints/basisset.h>
#include <libmints/onebody.h>
#include <libmints/twobody.h>
#include <libmints/integral.h>
#include <libmints/molecule.h>
#include <libscf_solver/integrator.h>

using namespace boost;
using namespace std;
using namespace psi;

namespace psi { namespace sapt {

void SAPT::lr_ints()
{
    //Open new LR Integrals file
    psio_->open(PSIF_SAPT_LRINTS,PSIO_OPEN_NEW);

    //Some parameters
    int norbs = calc_info_.nso;
    int naux = calc_info_.nri;
    int noccA = calc_info_.noccA;
    int noccB = calc_info_.noccB;
    int nvirA = calc_info_.nvirA;
    int nvirB = calc_info_.nvirB;
    double **CA = calc_info_.CA;
    double **CB = calc_info_.CB;

    //Prestripe
    psio_address lr_prestripe = PSIO_ZERO;
    double* bufferA = init_array(noccA*nvirA);
    for (int Q = 0; Q<naux; Q++)
        psio_->write(PSIF_SAPT_LRINTS,"A LR Integrals",(char *)&(bufferA[0]),sizeof(double)*noccA*nvirA,lr_prestripe,&lr_prestripe);
    free(bufferA);

    lr_prestripe = PSIO_ZERO;
    double* bufferB = init_array(noccB*nvirB);
    for (int Q = 0; Q<naux; Q++)
        psio_->write(PSIF_SAPT_LRINTS,"B LR Integrals",(char *)&(bufferB[0]),sizeof(double)*noccB*nvirB,lr_prestripe,&lr_prestripe);
    free(bufferB);

    //Setup Integrator
    
 
    //Setup Functional

    //Integrate and transform
    
    //Close LR Integrals file
    psio_->close(PSIF_SAPT_LRINTS,1);    
}

}}
