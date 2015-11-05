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

/*! \defgroup DETCI detci: The Determinant CI code */

/*! \file
    \ingroup DETCI
    \brief Determinant-based CI program

   DETCI

   DETERMINANT CI Program, incorporating Abelian point-group symmetry

   C. David Sherrill
   Center for Computational Quantum Chemistry
   University of Georgia
   August 1994

   Updated 3/95 to do frozen core and virtuals correctly
   Updated 5/95 to do RAS CI's again
   Updated 2/96 to clean up code and rename DETCI

*/

#include <cstdio>
#include <libmints/mints.h>
#include <pthread.h>
#include "structs.h"
#include "globals.h"
#include "tpool.h"
#include "civect.h"
#include "ciwave.h"

namespace psi { namespace detci {

PsiReturnType detci(Options &options);

}} // namespace psi::detci


namespace psi { namespace detci {


PsiReturnType detci(Options &options)
{

   boost::shared_ptr<Wavefunction> refwfn = Process::environment.wavefunction();
   boost::shared_ptr<CIWavefunction> ciwfn(new CIWavefunction(refwfn, options));

   if (Parameters.nthreads > 1)
     tpool_init(&thread_pool, Parameters.nthreads, CalcInfo.num_alp_str, 0);
                                /* initialize thread pool */

   if (Parameters.istop) {      /* Print size of space, other stuff, only   */
     ciwfn->cleanup();
     Process::environment.globals["CURRENT ENERGY"] = 0.0;
     Process::environment.globals["CURRENT CORRELATION ENERGY"] = 0.0;
     Process::environment.globals["CI TOTAL ENERGY"] = 0.0;
     Process::environment.globals["CI CORRELATION ENERGY"] = 0.0;

     return Success;
   }

   // MCSCF is special, we let it handle a lot of its own issues
   if (Parameters.mcscf){
     ciwfn->compute_mcscf();
   }
   else{
     // Transform and set ci integrals
     ciwfn->transform_ci_integrals();

     if (Parameters.mpn){
       ciwfn->compute_mpn();
       }
     else if (Parameters.cc)
       ciwfn->compute_cc();
     else
       ciwfn->diag_h();
   }

   // Finished CI, setting wavefunction parameters
   if(!Parameters.zaptn & Parameters.opdm){
     ciwfn->form_opdm();
     ciwfn->set_opdm();
   }
   else{
     ciwfn->set_opdm(true);
   }

   if (Parameters.tpdm) ciwfn->form_tpdm();
   if (Parameters.nthreads > 1) tpool_destroy(thread_pool, 1);
   if (Parameters.print_lvl > 0){
     outfile->Printf("\t\t \"A good bug is a dead bug\" \n\n");
     outfile->Printf("\t\t\t - Starship Troopers\n\n");
     outfile->Printf("\t\t \"I didn't write FORTRAN.  That's the problem.\"\n\n");
     outfile->Printf("\t\t\t - Edward Valeev\n\n");
   }


   ciwfn->cleanup();
   Process::environment.set_wavefunction((static_cast<boost::shared_ptr<Wavefunction> > (ciwfn)));
   return Success;
}



}} // namespace psi::detci


