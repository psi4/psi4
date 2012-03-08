#include<libplugin/plugin.h>
#include<libmints/matrix.h>
#include<libpsio/psio.hpp>
#include"cim.h"

INIT_PLUGIN

namespace psi{
  void RunCoupledCluster(Options &options,boost::shared_ptr<psi::CIM> wfn);
};


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
     options.add_double("THRESH2", 0.05);
     /*- cim threshold 3 -*/
     options.add_double("THRESH3", 1e-6);
  }
  return true;
}

extern "C" PsiReturnType
plugin_libcim(Options &options)
{  
  // cim class.  
  boost::shared_ptr<CIM>  cim(new CIM(options));

  /*
   *  loop over each cluster, {I}
   */
  cim->localClmo = SharedMatrix (new Matrix("Local Clmo",cim->nso,cim->nmo));
  cim->localFock = SharedMatrix (new Matrix("Local Fock",cim->nmo,cim->nmo));
  cim->Rii       = SharedMatrix (new Matrix("Rii",cim->ndoccact,cim->ndoccact));
  cim->centralfac = (double*)malloc(cim->ndoccact*sizeof(double));
  double ecim = 0.0;
  for (int i=0; i<cim->ndoccact; i++){
      if (cim->skip[i]) continue;
      cim->nfzc_     = cim->nfzc;
      cim->ndoccact_ = cim->ndoccact;
      cim->nvirt_    = cim->nvirt;
      cim->nfzv_     = cim->nfzv;

      /*
       * TODO
       * build virtual space for cluster, {I}
       */
      cim->VirtualSpaces(i);

      /*
       * TODO
       * make a wavefunction that has the right dimensions for this cluster
       * after this, local variables should be set: nfzc_, ndoccact_, nvirt_, nfzv_
       */
      cim->ClusterWavefunction(i);

      /*
       * TODO
       * canonicalize occupied and virtual spaces for cluster, {I}
       * this is where orbitals get reordered, too. also, count
       * how many times each central orbital occurs.
       */
      cim->QuasiCanonicalOrbitals(i);

      /*
       * TODO
       * transform integrals to local spaces for cluster, {I}
       */
      cim->TransformIntegrals(i,cim->localClmo);//boys->Clmo);

      /*
       * run correled calculation in cluster, {I}
       */
      RunCoupledCluster(options,cim);
      ecim += Process::environment.globals["CCSD CORRELATION ENERGY"];
  }
  printf("ecim %20.12lf\n",ecim);
  fprintf(outfile,"ecim %20.12lf\n",ecim);

  return  Success;
} // end plugin_libcim


}} // end namespace psi

