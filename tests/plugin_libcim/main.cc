#include<libplugin/plugin.h>
#include<libmints/matrix.h>
#include"cim.h"

INIT_PLUGIN

namespace psi{ namespace plugin_libcim{
extern "C" int 
read_options(std::string name, Options &options)
{
  if (name == "PLUGIN_LIBCIM"|| options.read_globals()) {
     /*- Convergence for the localization procedure-*/
     options.add_double("BOYS_CONVERGENCE", 1.0e-6);
     /*- Maximum number of localization iterations iterations -*/
     options.add_int("BOYS_MAXITER", 100);
     /*- Desired number of threads -*/
     options.add_int("NUM_THREADS", 1);
     /*- cim threshold 1 -*/
     options.add_double("THRESH1", 0.01);
     /*- cim threshold 2 -*/
     options.add_double("THRESH2", 0.03);
     /*- cim threshold 3 -*/
     options.add_double("THRESH3", 1e-3);
  }
  return true;
}

extern "C" PsiReturnType
plugin_libcim(Options &options)
{  
  // cim class
  boost::shared_ptr<CIM> cim (new CIM(options));

  return  Success;
} // end plugin_libcim


}} // end namespace psi

