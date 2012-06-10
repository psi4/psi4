#include <libplugin/plugin.h>
#include <libpsio/psio.hpp>
#include<libtrans/integraltransform.h>
#include<libtrans/mospace.h>

INIT_PLUGIN

namespace psi{ namespace plugin_mp2sort{

extern "C" int 
read_options(std::string name, Options &options)
{
  if (name == "PLUGIN_MP2SORT"|| options.read_globals()) {
  }
  return true;
}

extern "C" PsiReturnType
plugin_mp2sort(Options &options)
{  

  fflush(outfile);
  fprintf(outfile,"\n");
  fprintf(outfile, "        ********************************************************\n");
  fprintf(outfile, "        *                                                      *\n");
  fprintf(outfile, "        *         MP2 Integral Transformation and Sort         *\n");
  fprintf(outfile, "        *                                                      *\n");
  fprintf(outfile, "        ********************************************************\n");
  fprintf(outfile,"\n\n");
  fflush(outfile);

  fprintf(outfile,"        ==> Transform (OV|OV) integrals <==\n");
  fprintf(outfile,"\n");
  boost::shared_ptr<psi::Wavefunction> ref = Process::environment.reference_wavefunction();
  std::vector<boost::shared_ptr<MOSpace> > spaces;
  spaces.push_back(MOSpace::all);
  IntegralTransform ints(ref, spaces, IntegralTransform::Restricted,
           IntegralTransform::IWLAndDPD, IntegralTransform::QTOrder, IntegralTransform::None, true);
  ints.set_keep_dpd_so_ints(1);
  ints.set_keep_iwl_so_ints(1);
  ints.initialize();
  ints.transform_tei(MOSpace::all, MOSpace::all, MOSpace::all, MOSpace::all);
  fprintf(outfile,"\n");  

  //fprintf(outfile,"        ==> Sort (OV|OV) integrals <==\n");
  //sort(options);

  return  Success;
} // end plugin_mp2sort


}} // end namespace psi

