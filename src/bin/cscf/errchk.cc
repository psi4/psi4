/*! \file
    \ingroup CSCF
    \brief Enter brief description of file here 
*/
/* $Log$
 * Revision 1.2  2002/03/25 02:17:36  janssen
 * Get rid of tmpl.  Use new naming scheme for libipv1 includes.
 *
 * Revision 1.1.1.1  2000/02/04 22:52:30  evaleev
 * Started PSI 3 repository
 *
 * Revision 1.2  1999/08/17 19:04:14  evaleev
 * Changed the default symmetric orthogonalization to the canonical
 * orthogonalization. Now, if near-linear dependencies in the basis are found,
 * eigenvectors of the overlap matrix with eigenvalues less than 1E-6 will be
 * left out. This will lead to num_mo != num_so, i.e. SCF eigenvector is no
 * longer a square matrix. Had to rework some routines in libfile30, and add some.
 * The progrem prints out a warning if near-linear dependencies are found. TRANSQT
 * and a whole bunch of other codes has to be fixed to work with such basis sets.
 *
 * Revision 1.1.1.1  1999/04/12 16:59:26  evaleev
 * Added a version of CSCF that can work with CINTS.
 * -Ed
 * */

static char *rcsid = "$Id: errchk.cc 3840 2008-02-22 21:35:18Z evaleev $";

#include <cstdio>
#include <libipv1/ip_lib.h>

namespace psi { namespace cscf {

void errchk(int errcod, char* token)
{
  if (errcod) {
    fprintf(stderr,"ERROR: %s\n",ip_error_message(errcod));
    fprintf(stderr,"TOKEN: %s\n",token);
    }
  }

}} // namespace psi::cscf
