#include<libplugin/plugin.h>

#include "ccsd.h" 

INIT_PLUGIN

namespace psi{
  void RunCoupledCluster(Options &options);
};

namespace psi{ namespace plugin_ccsd_serial{
extern "C" int 
read_options(std::string name, Options &options)
{
  if (name == "PLUGIN_CCSD_SERIAL"|| options.read_globals()) {
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
      /*- compute triples by default */
      options.add_bool("COMPUTE_TRIPLES", true);
      /*- cutoff for occupation of MP2 NO orbitals in (T) -*/
      options.add_double("VIRTUAL_CUTOFF", 1.0e-6);
      /*- triples by default use the full virtual space -*/
      options.add_bool("TRIPLES_USE_NOS", false);
      /*- number of threads for triples, not set by default -*/
      options.add_int("NUM_THREADS", 1);
      /*- generate density-fitted integrals so we can skip
          transqt2() and OutOfCoreSort(). default false */
      options.add_bool("DF_INTEGRALS",false);
      /*- SCS MP2, default true -*/
      options.add_bool("SCS_MP2", true);
      /*- SCS CCSD, default true -*/
      options.add_bool("SCS_CCSD", true);
      /*- opposite-spin scaling factor -*/
      options.add_double("MP2_SCALE_OS",1.20);
      /*- same-spin scaling factor -*/
      options.add_double("MP2_SCALE_SS",1.0/3.0);
      /*- oppposite-spin scaling factor -*/
      options.add_double("CC_SCALE_OS", 1.27);
      /*- same-spin scaling factor -*/
      options.add_double("CC_SCALE_SS",1.13);
  }
  return true;
}

extern "C" PsiReturnType
plugin_ccsd_serial(Options &options)
{  
  RunCoupledCluster(options);
  return  Success;
} // end plugin_ccsd


}} // end namespace psi

