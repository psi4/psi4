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
      Process::environment.reference_wavefunction();
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
        
      for (int i=0; i<nmor; i++) {
        for (int j=0; j<nmor; j++) {
          Pap[i][j] = moinfo.opdm[i + mo_offset][j + mo_offset]; 
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

    oe->set_title("CC Density OPDM (Spin densities are not correct)");

    oe->compute();
}

}}
