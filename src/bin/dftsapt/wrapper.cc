/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
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
 *@END LICENSE
 */

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

#include "dftsapt.h"
#include "infsapt.h"

using namespace boost;

namespace psi { namespace dftsapt {

PsiReturnType dftsapt(boost::shared_ptr<Wavefunction> dimer, 
                      boost::shared_ptr<Wavefunction> mA, 
                      boost::shared_ptr<Wavefunction> mB)
{
    tstart();

    boost::shared_ptr<DFTSAPT> sapt = DFTSAPT::build(dimer,mA,mB);
    sapt->compute_energy();

    tstop();

    return Success;
}

PsiReturnType asapt(boost::shared_ptr<Wavefunction> dimer, 
                      boost::shared_ptr<Wavefunction> mA, 
                      boost::shared_ptr<Wavefunction> mB,
                      boost::shared_ptr<Wavefunction> eA,
                      boost::shared_ptr<Wavefunction> eB)
{
    tstart();

    boost::shared_ptr<ASAPT> sapt = ASAPT::build(dimer,mA,mB,eA,eB);
    sapt->compute_energy();

    tstop();

    return Success;
}

PsiReturnType infsapt(boost::shared_ptr<Wavefunction> dimer, 
                      boost::shared_ptr<Wavefunction> mA, 
                      boost::shared_ptr<Wavefunction> mB)
{
    tstart();

    std::vector<boost::shared_ptr<Wavefunction> > m;
    m.push_back(mA);
    m.push_back(mB);

    boost::shared_ptr<INFSAPT> sapt = INFSAPT::build(dimer,m);
    sapt->compute_energy();

    tstop();

    return Success;
}

}}
