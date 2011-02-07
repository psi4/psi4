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
#include <libiwl/iwl.hpp>
#include <libqt/qt.h>

#include <libmints/mints.h>
#include "uks.h"


using namespace std;
using namespace psi;

namespace psi { namespace scf {

UKS::UKS(Options& options, shared_ptr<PSIO> psio, shared_ptr<Chkpt> chkpt) : UHF(options, psio, chkpt)
{

}

UKS::UKS(Options& options, shared_ptr<PSIO> psio) : UHF(options, psio)
{

}

UKS::~UKS()
{

}

double UKS::compute_energy()
{
    fprintf(outfile,"  In UKS by Rob Parrish.\n");
    fprintf(outfile,"  Skeleton Functional, Method not complete.\n");
    return -1.0;
}

}}
