#include <libplugin/plugin.h>

INIT_PLUGIN

namespace psi{
  void RunCoupledCluster(Options &options);
};

namespace psi{ namespace plugin_df_ccsd{
extern "C" int 
read_options(std::string name, Options &options)
{
  if (name == "PLUGIN_DF_CCSD"|| options.read_globals()) {
     /*- Wavefunction type !expert -*/
     options.add_str("WFN", "CCSD");
     /*- Convergence for the CC amplitudes-*/
     options.add_double("R_CONVERGENCE", 1.0e-7);
     /*- Maximum number of CC iterations -*/
     options.add_int("MAXITER", 50);
     /*- Desired number of DIIS vectors -*/
     options.add_int("DIIS_MAX_VECS", 8);
     /*- For GPU code, cap the amount of memory registerred with the GPU -*/
     options.add_int("MAX_MAPPED_MEMORY", 1000);
     /*- Compute triples contribution? */
     options.add_bool("COMPUTE_TRIPLES", true);
     /*- Use MP2 NOs to truncate virtual space for (T)? -*/
     options.add_bool("NAT_ORBS", false);
     /*- Cutoff for occupation of MP2 NO orbitals in (T) -*/
     options.add_double("OCC_TOLERANCE", 1.0e-6);
     /*- Desired number of threads. This will override OMP_NUM_THREADS in (T) -*/
     options.add_int("NUM_THREADS", 1);
     /*- Do SCS-MP2? -*/
     options.add_bool("SCS_MP2", false);
     /*- Do SCS-CCSD? -*/
     options.add_bool("SCS_CCSD", false);
     /*- Opposite-spin scaling factor for SCS-MP2 -*/
     options.add_double("MP2_SCALE_OS",1.20);
     /*- Same-spin scaling factor for SCS-MP2 -*/
     options.add_double("MP2_SCALE_SS",1.0/3.0);
     /*- Oppposite-spin scaling factor for SCS-CCSD -*/
     options.add_double("CC_SCALE_OS", 1.27);
     /*- Same-spin scaling factor for SCS-CCSD -*/
     options.add_double("CC_SCALE_SS",1.13);
  }
  return true;
}

extern "C" PsiReturnType
plugin_df_ccsd(Options &options)
{  
  RunCoupledCluster(options);
  return  Success;
}


}} // end namespace psi

