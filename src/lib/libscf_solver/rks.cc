#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <algorithm>
#include <vector>
#include <utility>

#include <psifiles.h>
#include <physconst.h>
#include <libciomr/libciomr.h>
#include <libpsio/psio.hpp>
#include <libchkpt/chkpt.hpp>
#include <libipv1/ip_lib.h>
#include <libiwl/iwl.hpp>
#include <libqt/qt.h>

#include <libmints/matrix.h>
#include "rks.h"


using namespace std;
using namespace psi;

namespace psi { namespace scf {
    
RKS::RKS(Options& options, shared_ptr<PSIO> psio, shared_ptr<Chkpt> chkpt) : RHF(options, psio, chkpt)
{

}

RKS::~RKS()
{

}

double RKS::compute_energy()
{
    fprintf(outfile,"  In RKS by Rob Parrish.\n");
	fprintf(outfile,"  Skeleton Functional, Method not complete.\n");
    return -2.0;
}

}}
