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

/*! \file
    \ingroup CCDENSITY
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <libciomr/libciomr.h>
#include <libiwl/iwl.h>
#include <libchkpt/chkpt.h>
#include <libdpd/dpd.h>
#include <libqt/qt.h>
#include <psifiles.h>
#include <libmints/mints.h>
#include "MOInfo.h"
#include "Params.h"
#include "Frozen.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccdensity {
#include <physconst.h>

#define IOFF_MAX 32641
#define INDEX(i,j) ((i>j) ? (ioff[(i)]+(j)) : (ioff[(j)]+(i)))

void dipole(void)
{
    /* Call OEProp here for each root opdm */
    boost::shared_ptr<OEProp> oe(new OEProp());
    boost::shared_ptr<Wavefunction> wfn =
      Process::environment.wavefunction();
    boost::shared_ptr<Matrix> Ca = wfn->Ca();

    Dimension nmopi = wfn->nmopi();
    Dimension frzvpi = wfn->frzvpi();

    wfn->nalphapi().print();
    wfn->nbetapi().print();
    wfn->doccpi().print();
    wfn->soccpi().print();

    SharedMatrix Pa(new Matrix("P alpha", Ca->colspi(), Ca->colspi()));
    int mo_offset = 0;
    for (int h = 0; h < Ca->nirrep(); h++) {
      int nmo = nmopi[h]; 
      int nfv = frzvpi[h]; 
      int nmor = nmo - nfv;
      if (!nmo || !nmor) continue;
      double** Pap = Pa->pointer(h);
        
      // Loop over QT, convert to Pitzer
      for (int i=0; i<nmor; i++) {
        for (int j=0; j<nmor; j++) {
          int I = moinfo.pitzer2qt[i+mo_offset];
          int J = moinfo.pitzer2qt[j+mo_offset];
          Pap[i][j] = moinfo.opdm[I][J];
        } 
      }
      mo_offset += nmor; 
    }

    Pa->scale(0.5);
    //Pa->print();

    oe->set_Da_mo(Pa);
    if (!params.ref == 0) {
    	oe->set_Db_mo(Pa);
    }

    oe->add("DIPOLE");
    oe->add("QUADRUPOLE");
    oe->add("MULLIKEN_CHARGES");
    oe->add("NO_OCCUPATIONS");

    // TODO: This section needs work to duplicate the generic CC DIPOLE X, etc.
    //  into the exact method and/or root, like in DETCI
    psi::fprintf(outfile, "\nCC Density OPDM (Spin densities are not correct)\n\n");
    oe->set_title("CC");

    oe->compute();

//  Comments so that autodoc utility will find these PSI variables
//
//  Process::environment.globals["CC DIPOLE X"] =
//  Process::environment.globals["CC DIPOLE Y"] =
//  Process::environment.globals["CC DIPOLE Z"] =
//  Process::environment.globals["CC QUADRUPOLE XX"] =
//  Process::environment.globals["CC QUADRUPOLE XY"] =
//  Process::environment.globals["CC QUADRUPOLE XZ"] =
//  Process::environment.globals["CC QUADRUPOLE YY"] =
//  Process::environment.globals["CC QUADRUPOLE YZ"] =
//  Process::environment.globals["CC QUADRUPOLE ZZ"] =

}

}}
