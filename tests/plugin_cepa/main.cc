#include <libplugin/plugin.h>
#include <libmints/wavefunction.h>

#include "cepa.h" 

INIT_PLUGIN

// gpu ccsd class
namespace psi{
  void RunCoupledPair(Options &options,boost::shared_ptr<psi::Wavefunction> wfn);
};

namespace psi{ namespace plugin_cepa{
extern "C" int 
read_options(std::string name, Options &options)
{
  if (name == "PLUGIN_CEPA"|| options.read_globals()) {
      /*- The amount of information printed
             to the output file.  not used -*/
      options.add_int("PRINT", 1);
      /*- The amount of information printed
             to the output file. not used -*/
      options.add_int("DEBUG", 0);
      /*- default convergence of amplitudes-*/
      options.add_double("R_CONVERGENCE", 1.0e-7);
      /*- default maximum iterations -*/
      options.add_int("MAXITER", 100);
      /*- default memory available (mb) -*/
      options.add_int("MEMORY", 2000);
      /*- default number of DIIS iterations -*/
      options.add_int("DIIS_MAX_VECS", 8);
      /*- for GPU code, cap the amount of memory registerred with the GPU -*/
      options.add_int("MAX_MAPPED_MEMORY", 1000);
      /*- number of threads for triples, not set by default -*/
      options.add_int("NUM_THREADS", 1);
      /*- generate density-fitted integrals so we can skip
          transqt2() and OutOfCoreSort(). default false -*/
      options.add_bool("DF_INTEGRALS",false);
      /*- SCS MP2, default true -*/
      options.add_bool("SCS_MP2", true);
      /*- SCS CEPA, default true -*/
      options.add_bool("SCS_CEPA", true);
      /*- opposite-spin scaling factor -*/
      options.add_double("MP2_SCALE_OS",1.20);
      /*- same-spin scaling factor -*/
      options.add_double("MP2_SCALE_SS",1.0/3.0);
      /*- oppposite-spin scaling factor -*/
      options.add_double("CC_SCALE_OS", 1.27);
      /*- same-spin scaling factor -*/
      options.add_double("CC_SCALE_SS",1.13);
      /*- which cepa, default cepa(0) -*/
      options.add_str("CEPA_LEVEL","CEPA0");
      /*- compute dipole moment? default false-*/
      options.add_bool("DIPMOM",false);

  }
  return true;
}

extern "C" PsiReturnType
plugin_cepa(Options &options)
{  
  boost::shared_ptr<psi::Wavefunction> ref = Process::environment.reference_wavefunction();
  RunCoupledPair(options,ref);
  return  Success;
} // end plugin_cepa


}} // end namespace psi

