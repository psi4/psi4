#include "psi4-dec.h"
#include <libplugin/plugin.h>
#include <libparallel/parallel.h>
#include <liboptions/liboptions.h>
#include <libmints/mints.h>
#include <libpsio/psio.h>
#include <libchkpt/chkpt.hpp>
#include <libqt/qt.h>
#include <libciomr/libciomr.h>


#include "globals.h"
#include "ccsd.h" 
#include "runcoupledcluster.h"

INIT_PLUGIN

namespace psi{ namespace plugin_ccsd_serial_test{
extern "C" int 
read_options(std::string name, Options &options)
{
  if (name == "PLUGIN_CCSD_SERIAL_TEST"|| options.read_globals()) {
      /*- The amount of information printed
             to the output file -*/
      options.add_int("PRINT", 1);
      /*- The amount of information printed
             to the output file -*/
      options.add_int("DEBUG", 0);
      /*- default convergence of amplitudes-*/
      options.add_int("CONVERGENCE", 7);
      /*- default maximum iterations -*/
      options.add_int("MAXITER", 100);
      /*- default memory available (mb) -*/
      options.add_int("MEMORY", 2000);
      /*- default number of DIIS iterations */
      options.add_int("MAX_DIIS_VECS", 8);
      /*- default number of DIIS iterations */
      options.add_int("MAX_MAPPED_MEMORY", 1000);
      /*- do not compute triples by default */
      options.add_bool("COMPUTE_TRIPLES", false);
  }
  return true;
}

extern "C" PsiReturnType
plugin_ccsd_serial_test(Options &options)
{  
  RunCoupledCluster(options);
  return  Success;
} // end plugin_ccsd


}} // end namespace psi

