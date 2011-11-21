#include"psi4-dec.h"
#include"ccsd.h"
#include <libplugin/plugin.h>
#include<libdpd/dpd.h>
#include<boost/shared_ptr.hpp>
#include<liboptions/liboptions.h>
#include <ccfiles.h>
#include <libpsio/psio.hpp>
#include <libciomr/libciomr.h>
#include <ccfiles.h>

using namespace psi;
using namespace boost;

// Forward declaration to call cctriples
namespace psi { namespace cctriples {
PsiReturnType cctriples(Options &options);
}}
// Forward declaration to call ccsort
namespace psi { namespace ccsort {
PsiReturnType ccsort(Options &options);
}}
// Forward declaration to call ccenergy
namespace psi { namespace ccenergy {
  PsiReturnType ccenergy(Options &options);
  struct dpd_file4_cache_entry *priority_list(void);
  int **cacheprep_rhf(int level, int *cachefiles);
}}


namespace psi{
PsiReturnType triples(boost::shared_ptr<psi::CoupledCluster>ccsd,Options&options);

PsiReturnType triples(boost::shared_ptr<psi::CoupledCluster>ccsd,Options&options){

  // call ccsort so cctriples will have the integrals it needs
  psi::ccsort::ccsort(options);

  // grab parameters from CC_INFO
  boost::shared_ptr<PSIO> psio(new PSIO());
  psio->open(CC_INFO,PSIO_OPEN_OLD);
  int nactive;
  psio->read_entry(CC_INFO, "No. of Active Orbitals", (char *) &(nactive),sizeof(int));
  int*occ_sym = init_int_array(nactive);
  int*vir_sym = init_int_array(nactive);
  psio->read_entry(CC_INFO, "Active Occ Orb Symmetry",
          (char *) occ_sym, sizeof(int)*nactive);
  psio->read_entry(CC_INFO, "Active Virt Orb Symmetry",
          (char *) vir_sym, sizeof(int)*nactive);
  int*occpi = init_int_array(ccsd->nirreps);
  int*virtpi = init_int_array(ccsd->nirreps);
  psio->read_entry(CC_INFO, "Active Occ Orbs Per Irrep",
          (char *) occpi, sizeof(int)*ccsd->nirreps);
  psio->read_entry(CC_INFO, "Active Virt Orbs Per Irrep",
          (char *) virtpi, sizeof(int)*ccsd->nirreps);
  psio->close(CC_INFO,1);
  psio.reset();

  psio_open(CC_INFO,PSIO_OPEN_OLD);
  psio_open(CC_OEI,PSIO_OPEN_OLD);
  psio_open(CC_TAMPS,PSIO_OPEN_OLD);

  // write ccsd energy to disk - cctriples won't run without it
  psio_write_entry(CC_INFO, "CCSD Energy", (char *) &(ccsd->energy),sizeof(double));

  // set up cache stuff like ccenergy would so we can use dpd the
  // same exact way as it's used in there.
  int **cachelist, *cachefiles;
  struct dpd_file4_cache_entry *priority;
  cachefiles = init_int_array(PSIO_MAXUNIT);
  int cachelev = options.get_int("CACHELEV");
  std::string tempcachetype = "";
  int cachetype = 1;
  tempcachetype = options.get_str("CACHETYPE");
  if(tempcachetype == "LOW") cachetype = 1;
  else if(tempcachetype == "LRU") cachetype = 0;
  else
    throw PsiException("Error in input: invalid CACHETYPE", __FILE__, __LINE__);

  cachelist = psi::ccenergy::cacheprep_rhf(cachelev, cachefiles);
  priority  = psi::ccenergy::priority_list();

  long int memory = Process::environment.get_memory();

  // dpd!
  dpd_init(0, ccsd->nirreps, memory, cachetype, cachefiles,
       cachelist, priority, 2, occpi, occ_sym, virtpi, vir_sym);

  long int o = ccsd->ndoccact;
  long int v = ccsd->nvirt;

  // t1 amps
  dpdfile2 t1;
  dpd_file2_init(&t1, CC_OEI, 0, 0, 1, "tIA");
  dpd_file2_mat_init(&t1);
  dpd_file2_mat_rd(&t1);

  for (long int i=0; i<t1.params->rowtot[0]; i++){
      for (long int a=0; a<t1.params->coltot[0]; a++){
          t1.matrix[0][i][a] = ccsd->t1[a*o+i];
      }
  }
  dpd_file2_mat_wrt(&t1);
  dpd_file2_mat_close(&t1);
  dpd_file2_close(&t1);

  // t2 amps
  dpdbuf4 t2;
  dpd_buf4_init(&t2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
  dpd_buf4_mat_irrep_init(&t2, 0);
  dpd_buf4_mat_irrep_rd(&t2, 0);
  for (long int ij=0; ij< t2.params->rowtot[0]; ij++){
      long int i = t2.params->roworb[0][ij][0];
      long int j = t2.params->roworb[0][ij][1];
      for (long int ab=0; ab< t2.params->coltot[0]; ab++){
          long int a = t2.params->colorb[0][ab][0];
          long int b = t2.params->colorb[0][ab][1];
          t2.matrix[0][ij][ab] = ccsd->tb[a*o*o*v+b*o*o+i*o+j];

      }
  }
  dpd_buf4_mat_irrep_wrt(&t2, 0);
  dpd_buf4_mat_irrep_close(&t2, 0);
  dpd_buf4_close(&t2);

  dpd_close(0);

  psio_close(CC_INFO,1);
  psio_close(CC_OEI,1);
  psio_close(CC_TAMPS,1);

  // finally call crawford's triples code
  return psi::cctriples::cctriples(options);

}

} // end of namespace



