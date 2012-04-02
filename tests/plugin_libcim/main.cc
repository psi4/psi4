#include<libplugin/plugin.h>
#include<libmints/matrix.h>
#include<libpsio/psio.hpp>
#include<lib3index/dftensor.h>
#include<libciomr/libciomr.h>
#include"cim.h"
#include"blas.h"

INIT_PLUGIN

namespace psi{
  void RunCoupledCluster(Options &options,boost::shared_ptr<psi::CIM> wfn);
  void RunCoupledPair(Options &options,boost::shared_ptr<psi::CIM> wfn);
};

namespace psi{ namespace plugin_libcim{
extern "C" int 
read_options(std::string name, Options &options)
{
  if (name == "PLUGIN_LIBCIM"|| options.read_globals()) {
     /*- Convergence for the localization procedure -*/
     options.add_double("BOYS_CONVERGENCE", 1.0e-6);
     /*- Maximum number of iterations to converge the orbital localization
     procedure -*/
     options.add_int("BOYS_MAXITER", 100);
     /*- Desired number of threads -*/
     options.add_int("NUM_THREADS", 1);
     /*- decim tolerance 1 -*/
     options.add_double("CIM_DE_TOLERANCE1", 0.01);
     /*- decim tolerance 2 -*/
     options.add_double("CIM_DE_TOLERANCE2", 0.05);
     /*- secim tolerance -*/
     options.add_double("CIM_SE_TOLERANCE", 0.001);
     /*- Minimum occupation (eigenvalues of the MP2 OPDM) 
         below which virtual natural orbitals are discarded 
         for for a given cluster -*/
     options.add_double("OCC_TOLERANCE", 5e-5);
     /*- Should the CIM procedure return after the occupied domains are
     determined?  This parameter is used internally by the python driver
     if the calculation is going to be run in parallel.  Changing this
     won't have any effect on the procedure. -*/
     options.add_bool("CIM_INITIALIZE", false);
     /*- For which cluster number should the correlation energy be
     evaluated? This parameter is used internally by the python driver if
     the calculation is going to be run in parallel.  Changing this won't
     have any effect on the procedure. -*/
     options.add_int("CIM_CLUSTER_NUM", 0);
     /*- CIM central domain type (dual- or single-environment CIM)
     Recommeneded SECIM for all CIM computations. -*/
     options.add_str("CIM_DOMAIN_TYPE", "SECIM", "SECIM DECIM");
     /*- Maximum error allowed (Max error norm in Delta tensor) in the
     approximate energy denominators employed in evaluating the scaled
     opposite-spin MP2 OPDM used in defining the virtual orbitals for each
     CIM cluster.  The default may be more conservative than is necessary
     in practice. -*/
     options.add_double("DENOMINATOR_DELTA", 1.0E-6);
    /*- The correlated to be used in the CIM procedure. This parameter 
        is used internally by the python driver. Changing its value 
        won’t have any effect on the procedure -*/
    options.add_str("CIM_CORRELATED_METHOD", "CCSD(T)");
  }
  return true;
}

extern "C" PsiReturnType
plugin_libcim(Options &options)
{  
  // cim class.  
  boost::shared_ptr<CIM>  cim(new CIM());

  fflush(outfile);
  fprintf(outfile,"\n\n");
  fprintf(outfile, "        *******************************************************\n");
  fprintf(outfile, "        *                                                     *\n");
  fprintf(outfile, "        *                 Cluster-in-molecule                 *\n");
  fprintf(outfile, "        *                   Eugene DePrince                   *\n");
  fprintf(outfile, "        *                                                     *\n");
  fprintf(outfile, "        *******************************************************\n");
  fprintf(outfile,"\n\n");
  fflush(outfile);

  cim->options_ = Process::environment.options;
  cim->BuildClusters();

  if (options.get_bool("CIM_INITIALIZE")){
     Process::environment.globals["CIM CLUSTERS"] = cim->ndomains;
     return Success;
  }

  /*
   *  3-index integrals to approximate mp2 density:
   */
  int nocc = cim->doccpi()[0];
  int nvir = cim->nmopi()[0]-cim->doccpi()[0];
  int aocc = nocc-cim->frzcpi()[0];
  int avir = nvir-cim->frzvpi()[0];
  // build 3-index tensors (in lmo basis)
  boost::shared_ptr<BasisSetParser> parser(new Gaussian94BasisSetParser());
  boost::shared_ptr<BasisSet> auxilliary = BasisSet::construct(parser, cim->molecule(), "DF_BASIS_MP2");
  boost::shared_ptr<DFTensor> DF (new DFTensor(cim->basisset(),auxilliary,cim->boys->Clmo,nocc,nvir,aocc,avir,options));
  SharedMatrix tmpQov = DF->Qov();
  cim->Qov = tmpQov->pointer();
  cim->nQ = auxilliary->nbf();

  /*
   *  loop over each cluster, {I}
   */
  cim->localClmo = SharedMatrix (new Matrix("Local Clmo",cim->nso,cim->nmo));
  cim->localFock = SharedMatrix (new Matrix("Local Fock",cim->nmo,cim->nmo));
  cim->Rii       = SharedMatrix (new Matrix("Rii",cim->ndoccact,cim->ndoccact));
  cim->centralfac = (double*)malloc(cim->ndoccact*sizeof(double));
  double ecim = 0.0;

  int clusternum = -1;
  for (int i=0; i<cim->maxndomains; i++){
      if (cim->skip[i]) continue;

      clusternum++;

      if (options.get_int("CIM_CLUSTER_NUM")!=clusternum) continue;

      cim->nfzc_     = cim->nfzc;
      cim->ndoccact_ = cim->ndoccact;
      cim->nvirt_    = cim->nvirt;
      cim->nfzv_     = cim->nfzv;

      fprintf(outfile,"\n");
      fprintf(outfile,"  ==> Cluster {%i} <==\n",clusternum);
      fprintf(outfile,"\n");
      fprintf(outfile,"        Number of occupied orbitals in original space:  %5i\n",cim->ndoccact_);
      fprintf(outfile,"        Number of occupied orbitals in truncated space: %5i\n",cim->ncentral[i]+cim->nmodomain[i]+cim->nenv[i]);
      fprintf(outfile,"\n");

      /*
       * build virtual space for cluster, {I}
       */
      tstart();
      cim->VirtualSpaces(i,clusternum);
      tstop();

      /*
       * make a wavefunction that has the right dimensions for this cluster
       * after this, local variables should be set: 
       * nfzc_, ndoccact_, nvirt_, nfzv_
       */
      cim->ClusterWavefunction(i);

      /*
       * canonicalize occupied and virtual spaces for cluster, {I}
       * this is where orbitals get reordered, too. also, count
       * how many times each central orbital occurs.
       */
      cim->QuasiCanonicalOrbitals(i);

      /*
       * transform integrals to local spaces for cluster, {I}
       */
      tstart();
      cim->TransformIntegrals(i,cim->localClmo,clusternum);
      tstop();

      /*
       * run correled calculation in cluster, {I}
       */
      if (options.get_str("CIM_CORRELATED_METHOD") == "CCSD" || options.get_str("CIM_CORRELATED_METHOD") == "CCSD(T)"){
         RunCoupledCluster(options,cim);
         double dum = Process::environment.globals["CCSD CORRELATION ENERGY"];
         Process::environment.globals["CURRENT CLUSTER CCSD CORRELATION ENERGY"] = dum;
         
         fprintf(outfile,"\n");
         fprintf(outfile,"        Cluster {%i} CCSD energy contribution:     %20.12lf\n",clusternum,dum);
         if (options.get_bool("COMPUTE_TRIPLES")){
            double dumt = Process::environment.globals["(T) CORRELATION ENERGY"];
            Process::environment.globals["CURRENT CLUSTER (T) CORRELATION ENERGY"] = dumt;
            Process::environment.globals["CURRENT CLUSTER CCSD(T) CORRELATION ENERGY"] = dum + dumt;
            fprintf(outfile,"        Cluster {%i} (T) energy contribution:      %20.12lf\n",clusternum,dumt);
            fprintf(outfile,"        Cluster {%i} CCSD(T) energy contribution:  %20.12lf\n",clusternum,dum+dumt);
         }
         fprintf(outfile,"\n");
         fflush(outfile);
      }
      else if (options.get_str("CIM_CORRELATED_METHOD") == "CEPA"){
         RunCoupledPair(options,cim);
         double dum = Process::environment.globals["CURRENT CORRELATION ENERGY"];
         Process::environment.globals["CURRENT CLUSTER CEPA CORRELATION ENERGY"] = dum;
         
         fprintf(outfile,"\n");
         fprintf(outfile,"        Cluster {%i} correlation energy contribution:     %20.12lf\n",clusternum,dum);
         fprintf(outfile,"\n");
         fflush(outfile);
      }

      /*
       *  need to replace eps in wavefunction. this is stupid.
       */
      double*eps = cim->epsilon_a()->pointer();
      F_DCOPY(cim->nmo,cim->epsSave,1,eps,1);
      free(cim->epsSave);

      return Success;
  }
  /*double escf    = Process::environment.globals["SCF TOTAL ENERGY"];
  fprintf(outfile,"\n");
  fprintf(outfile,"        CIM correlation energy: %20.12lf\n",ecim);
  fprintf(outfile,"      * CIM total energy:       %20.12lf\n",ecim+escf);
  fprintf(outfile,"\n");
  fflush(outfile);*/

  return  Success;
} // end plugin_libcim


}} // end namespace psi

