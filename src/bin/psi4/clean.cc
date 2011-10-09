/*! \file clean
    \defgroup PSI4
*/

#include <boost/shared_ptr.hpp>
#include <libpsio/psio.h>
#include <libpsio/psio.hpp>

namespace psi {

    extern FILE *outfile;

/*!
** psiclean():
**
** Remove all psio files in the event of a crash.
** Uses the "psi.clean" file
*/

void psiclean(void) {

PSIOManager::shared_object()->crashclean();

}

} //end ::psi

