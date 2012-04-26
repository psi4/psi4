#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <sstream>

#include <psifiles.h>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>
#include <libchkpt/chkpt.h>
#include <libpsio/psio.hpp>
#include <libchkpt/chkpt.hpp>
#include <libiwl/iwl.h>
#include <libqt/qt.h>
#include <libmints/mints.h>
#include <psi4-dec.h>

#include "mp2.h"

using namespace boost;

namespace psi { namespace dfmp2 {

PsiReturnType dfmp2(Options & options)
{
    tstart();

    boost::shared_ptr<PSIO> psio(new PSIO);
    boost::shared_ptr<Chkpt> chkpt(new Chkpt(psio, PSIO_OPEN_OLD));

    boost::shared_ptr<DFMP2> dfmp2;
    if (options.get_str("REFERENCE") == "RHF" || options.get_str("REFERENCE") == "RKS") {
        dfmp2 = boost::shared_ptr<DFMP2>(new RDFMP2(options,psio,chkpt)); 
    } else if (options.get_str("REFERENCE") == "UHF" || options.get_str("REFERENCE") == "UKS") {
        dfmp2 = boost::shared_ptr<DFMP2>(new UDFMP2(options,psio,chkpt)); 
    } else if (options.get_str("REFERENCE") == "ROHF") {
        dfmp2 = boost::shared_ptr<DFMP2>(new RODFMP2(options,psio,chkpt)); 
    } else {
        throw PSIEXCEPTION("DFMP2: Unrecognized reference");
    }
    dfmp2->compute_energy();

    tstop();

    return Success;
}

PsiReturnType dfmp2grad(Options & options)
{
    tstart();

    boost::shared_ptr<PSIO> psio(new PSIO);
    boost::shared_ptr<Chkpt> chkpt(new Chkpt(psio, PSIO_OPEN_OLD));

    boost::shared_ptr<DFMP2> dfmp2;
    if (options.get_str("REFERENCE") == "RHF" || options.get_str("REFERENCE") == "RKS") {
        dfmp2 = boost::shared_ptr<DFMP2>(new RDFMP2(options,psio,chkpt)); 
    } else if (options.get_str("REFERENCE") == "UHF" || options.get_str("REFERENCE") == "UKS") {
        dfmp2 = boost::shared_ptr<DFMP2>(new UDFMP2(options,psio,chkpt)); 
    } else if (options.get_str("REFERENCE") == "ROHF") {
        dfmp2 = boost::shared_ptr<DFMP2>(new RODFMP2(options,psio,chkpt)); 
    } else {
        throw PSIEXCEPTION("DFMP2: Unrecognized reference");
    }
    dfmp2->compute_gradient();

    tstop();

    return Success;
}

}}
