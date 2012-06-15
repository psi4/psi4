#include"psi4-dec.h"
#include <libpsio/psio.hpp>
#include<libtrans/integraltransform.h>
#include<libtrans/mospace.h>
#include<libmints/wavefunction.h>
#include<libmints/matrix.h>

#ifdef _OPENMP
    #include<omp.h>
#else
    #define omp_get_wtime() 0.0
#endif

#include"coupledpair.h"
#include"blas.h"

using namespace psi;
namespace psi{ namespace cepa{
  long int Position(long int i,long int j);
}}

namespace psi{ namespace cepa{

void TransformOXXX(boost::shared_ptr<Wavefunction>wfn);
void TransformAllIntegrals(boost::shared_ptr<Wavefunction>wfn);

void TransformIntegrals(boost::shared_ptr<Wavefunction>wfn,Options&options){
  if (options.get_bool("CEPA_VABCD_DIRECT")){
     TransformOXXX(wfn);
  }else{
     TransformAllIntegrals(wfn);
  }
}

void TransformAllIntegrals(boost::shared_ptr<Wavefunction>wfn){
  fflush(outfile);
  fprintf(outfile,"\n");
  fprintf(outfile, "        *******************************************************\n");
  fprintf(outfile, "        *                                                     *\n");
  fprintf(outfile, "        *               Integral Transformation               *\n");
  fprintf(outfile, "        *                                                     *\n");
  fprintf(outfile, "        *******************************************************\n");
  fprintf(outfile,"\n\n");
  fflush(outfile);

  boost::shared_ptr<PSIO> psio(_default_psio_lib_);
  std::vector<boost::shared_ptr<MOSpace> > spaces;

  spaces.push_back(MOSpace::all);
  IntegralTransform ints(wfn, spaces, IntegralTransform::Restricted,
           IntegralTransform::IWLOnly, IntegralTransform::QTOrder,
           IntegralTransform::None, false);
  ints.set_keep_dpd_so_ints(1);
  ints.set_keep_iwl_so_ints(1);
  ints.initialize();

  fprintf(outfile,"        ==> Transform (xx|xx) integrals <==\n");
  fprintf(outfile,"\n");
  ints.transform_tei(MOSpace::all, MOSpace::all, MOSpace::all, MOSpace::all);
  fprintf(outfile,"\n");
}
void TransformOXXX(boost::shared_ptr<Wavefunction>wfn){
  fflush(outfile);
  fprintf(outfile,"\n");
  fprintf(outfile, "        *******************************************************\n");
  fprintf(outfile, "        *                                                     *\n");
  fprintf(outfile, "        *               Integral Transformation               *\n");
  fprintf(outfile, "        *                                                     *\n");
  fprintf(outfile, "        *******************************************************\n");
  fprintf(outfile,"\n\n");
  fflush(outfile);

  boost::shared_ptr<PSIO> psio(_default_psio_lib_);
  std::vector<boost::shared_ptr<MOSpace> > spaces;

  spaces.push_back(MOSpace::occ);
  spaces.push_back(MOSpace::all);
  IntegralTransform ints(wfn, spaces, IntegralTransform::Restricted,
           IntegralTransform::IWLOnly, IntegralTransform::QTOrder,
           IntegralTransform::None, false);
  ints.set_keep_dpd_so_ints(1);
  ints.set_keep_iwl_so_ints(1);
  ints.initialize();

  fprintf(outfile,"        ==> Transform (ox|xx) integrals <==\n");
  fprintf(outfile,"\n");
  ints.transform_tei(MOSpace::occ, MOSpace::all, MOSpace::all, MOSpace::all);
  fprintf(outfile,"\n");
}

}} // end of namespace
