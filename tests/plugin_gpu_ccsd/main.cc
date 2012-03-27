#include<libplugin/plugin.h>
#include<libciomr/libciomr.h>
#include "gpu_ccsd.h" 

#define HAVE_GPU

INIT_PLUGIN

namespace psi{
  PsiReturnType triples(boost::shared_ptr<psi::GPUCoupledCluster>ccsd,Options&options);
}

namespace psi{ namespace plugin_gpu_ccsd{
extern "C" int 
read_options(std::string name, Options &options)
{
  if (name == "PLUGIN_GPU_CCSD"|| options.read_globals()) {
     /*- Memory in mb.  Use this to override global option, 
         which is limited to <16gb because of dpd in transqt2 -*/
     options.add_int("MEMORY", 256);
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
plugin_gpu_ccsd(Options &options)
{  
  tstart();

  // pointer to gpu ccsd class
  boost::shared_ptr<GPUCoupledCluster> ccsd(new GPUCoupledCluster);

  ccsd->WriteBanner(options);
  ccsd->Initialize(options);
  ccsd->AllocateMemory();
  ccsd->ReadIntegrals();
  PsiReturnType status = ccsd->CCSDIterations(options);

  // mp2 energy
  Process::environment.globals["MP2 CORRELATION ENERGY"] = ccsd->emp2;
  Process::environment.globals["MP2 TOTAL ENERGY"] = ccsd->emp2 + ccsd->escf;

  // ccsd energy
  Process::environment.globals["CCSD CORRELATION ENERGY"] = ccsd->eccsd;
  Process::environment.globals["CCSD TOTAL ENERGY"] = ccsd->eccsd + ccsd->escf;
  Process::environment.globals["CURRENT ENERGY"] = ccsd->eccsd + ccsd->escf;

  tstop();

  if (options.get_bool("COMPUTE_TRIPLES")){
     tstart();

     // triples
     status = psi::triples(ccsd,options);
     if (status == Failure){
        throw PsiException(
           "Whoops, the (T) correction died.",__FILE__,__LINE__);
     }

     // ccsd(t) energy
     Process::environment.globals["(T) CORRELATION ENERGY"] = ccsd->et;
     Process::environment.globals["CCSD(T) CORRELATION ENERGY"] = ccsd->eccsd + ccsd->et;
     Process::environment.globals["CCSD(T) TOTAL ENERGY"] = ccsd->eccsd + ccsd->et + ccsd->escf;
     Process::environment.globals["CURRENT ENERGY"] = ccsd->eccsd + ccsd->et + ccsd->escf;
     tstop();
  }

  ccsd->CudaFinalize();

  return  Success;
} // end plugin_gpu_ccsd


}} // end namespace psi

