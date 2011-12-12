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
#include "dfmp2.h"
//#include "mad_mp2.h"

#include <psi4-dec.h>

using namespace boost;

namespace psi { namespace dfmp2 {

std::string to_string(const int val);   // In matrix.cpp

PsiReturnType dfmp2(Options & options)
{
    tstart();

    boost::shared_ptr<PSIO> psio(new PSIO);


#if HAVE_MADNESS
    mad_mp2::MAD_MP2 madmp2(options, psio);
    madmp2.compute_energy();
#else
    boost::shared_ptr<Chkpt> chkpt(new Chkpt(psio, PSIO_OPEN_OLD));
    DFMP2 df(options, psio, chkpt);
    df.compute_energy();
#endif
    // Shut down psi.

    tstop();

    return Success;
}

}}
