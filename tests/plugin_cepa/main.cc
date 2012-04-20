#include <libplugin/plugin.h>
#include <libmints/wavefunction.h>

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
      /*- Desired convergence for the t1 and t2 amplitudes, defined as
      the norm of the change in the amplitudes between iterations.-*/
      options.add_double("R_CONVERGENCE", 1.0e-7);
      /*- Maximum number of iterations to converge the t1 and t2
      amplitudes. -*/
      options.add_int("MAXITER", 100);
      /*- Number of vectors to store for DIIS extrapolation. -*/
      options.add_int("DIIS_MAX_VECS", 8);
      /*- Opposite-spin scaling factor for SCS-MP2. -*/
      options.add_double("MP2_SCALE_OS",1.20);
      /*- Same-spin scaling factor for SCS-MP2-*/
      options.add_double("MP2_SCALE_SS",1.0/3.0);
      /*- Perform SCS-CEPA? If true, note that the
      default values for the spin component scaling factors
      are optimized for the CCSD method. -*/
      options.add_bool("SCS_CEPA", false);
      /*- Oppposite-spin scaling factor for SCS-CEPA. -*/
      options.add_double("CEPA_SCALE_OS", 1.27);
      /*- Same-spin scaling factor for SCS-CEPA. -*/
      options.add_double("CEPA_SCALE_SS",1.13);
      /*- Which coupled-pair method is called?  This parameter is
      used internally by the python driver.  Changing its value 
      won't have any effect on the procedure. -*/
      options.add_str("CEPA_LEVEL","CEPA0");
      /*- Compute the dipole moment? Note that quadrapole moments
      will also be computed if PRINT >= 2. -*/
      options.add_bool("DIPMOM",false);
      /*- Flag to exclude singly excited configurations from the
      computation. Note that this algorithm is not optimized for
      doubles-only computations. -*/
      options.add_bool("CEPA_NO_SINGLES",false);
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

