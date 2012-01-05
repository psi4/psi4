#include "integraltransform.h"
#include <libchkpt/chkpt.hpp>
#include <libpsio/psio.hpp>
#include <libciomr/libciomr.h>
#include <libiwl/iwl.hpp>
#include <libqt/qt.h>
#include <math.h>
#include <ctype.h>
#include <stdio.h>
#include "psifiles.h"
#include "ccfiles.h"
#include "mospace.h"
#define EXTERN
#include <libdpd/dpd.gbl>

using namespace boost;
using namespace psi;

/**
 * Transform the two-electron integrals from the SO to the MO basis in the spaces specified
 *
 * @param s1 - the MO space for the first index
 * @param s2 - the MO space for the second index
 * @param s3 - the MO space for the third index
 * @param s4 - the MO space for the fourth index
 */
void
IntegralTransform::transform_tei(const shared_ptr<MOSpace> s1, const shared_ptr<MOSpace> s2,
                                 const shared_ptr<MOSpace> s3, const shared_ptr<MOSpace> s4)
{
    check_initialized();

    /* The only difficulty here is that we have to figure out which integrals are unique,
     * which requires knowing the allowed permutations.  This is easy - it's just a bunch of
     * if statements, checking things like e.g. "if s1==s2" to see if bra indices can be permuted.
     * I haven't done this yet because I don't anticipate much demand for IWL output.  If
     * somebody need this, they can contact me for help to implement it if needed.
     * Andy Simmonett, 11/09*/
    if(useIWL_ && !(s1==s2 && s1==s3 && s1==s4))
        throw FeatureNotImplemented("libtrans", "mixed spaces with IWL output", __FILE__, __LINE__);

    transform_tei_first_half(s1, s2);
    transform_tei_second_half(s1, s2, s3, s4);
}
