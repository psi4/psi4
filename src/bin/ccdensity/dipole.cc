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
    boost::shared_ptr<Wavefunction> wfn = Process::environment.wavefunction();
    boost::shared_ptr<Matrix> Ca = wfn->Ca();
    boost::shared_ptr<Matrix> Cb = wfn->Cb();

//    Ca->print();

    Dimension nmopi = wfn->nmopi();
    Dimension frzvpi = wfn->frzvpi();

    wfn->nalphapi().print();
    wfn->nbetapi().print();
    wfn->doccpi().print();
    wfn->soccpi().print();
    outfile->Printf("Wfn name = %s\n", wfn->name().c_str());
    outfile->Printf("Same alpha/beta density? %d\n", wfn->same_a_b_dens());

//    outfile->Printf("Alpha OPDM:\n");
//    mat_print(moinfo.opdm_a, moinfo.nmo, moinfo.nmo, "outfile");
//    outfile->Printf("Beta OPDM:\n");
//    mat_print(moinfo.opdm_b, moinfo.nmo, moinfo.nmo, "outfile");

    SharedMatrix Pa(new Matrix("P alpha", Ca->colspi(), Ca->colspi()));
    SharedMatrix Pb(new Matrix("P beta", Cb->colspi(), Cb->colspi()));
    int mo_offset = 0;
    for (int h = 0; h < Ca->nirrep(); h++) {
      int nmo = nmopi[h];
      int nfv = frzvpi[h];
      int nmor = nmo - nfv;
      if (!nmo || !nmor) continue;

      // Loop over QT, convert to Pitzer
      double **Pap = Pa->pointer(h);
      double **Pbp = Pb->pointer(h);
      for (int i=0; i<nmor; i++) {
        for (int j=0; j<nmor; j++) {
          int I = moinfo.pitzer2qt[i+mo_offset];
          int J = moinfo.pitzer2qt[j+mo_offset];
          if(wfn->same_a_b_dens())
            Pap[i][j] = moinfo.opdm[I][J];
          else {
            Pap[i][j] = moinfo.opdm_a[I][J];
            Pbp[i][j] = moinfo.opdm_b[I][J];
          }
        }
      }
      mo_offset += nmo;
    }

//    Pa->print();
//    Pb->print();

    if(!wfn->same_a_b_dens()) {
      SharedMatrix Pa_V(new Matrix("Alpha Density Eigenvectors", Ca->colspi(), Ca->colspi()));
      SharedVector Pa_v(new Vector("Alpha Density Eigenvalues", Ca->colspi()));
      Pa->diagonalize(Pa_V, Pa_v, descending);
      Pa_v->print();
      SharedMatrix Pb_V(new Matrix("Beta Density Eigenvectors", Cb->colspi(), Cb->colspi()));
      SharedVector Pb_v(new Vector("Beta Density Eigenvalues", Cb->colspi()));
      Pb->diagonalize(Pb_V, Pb_v, descending);
      Pb_v->print();
    }

    if(wfn->same_a_b_dens()) Pa->scale(0.5);
    //Pa->print();

    oe->set_Da_mo(Pa);
    if(!wfn->same_a_b_dens()) oe->set_Db_mo(Pb);

    oe->add("DIPOLE");
    oe->add("QUADRUPOLE");
    oe->add("MULLIKEN_CHARGES");
    oe->add("NO_OCCUPATIONS");

    // TODO: This section needs work to duplicate the generic CC DIPOLE X, etc.
    //  into the exact method and/or root, like in DETCI
    outfile->Printf( "\nCC Density OPDM\n\n");
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
