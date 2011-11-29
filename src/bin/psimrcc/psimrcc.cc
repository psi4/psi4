/**
 *  @defgroup PSIMRCC PSIMRCC is a code for SR/MRCC computations
 *  @file psimrcc.cpp
 *  @ingroup (PSIMRCC)
 *  @brief Contains main() and global variables
*/

#include <liboptions/liboptions.h>
#include <libmoinfo/libmoinfo.h>
#include <libutil/libutil.h>

#include "blas.h"
#include "sort.h"
#include "mrcc.h"
#include "idmrpt2.h"
#include "mp2_ccsd.h"
#include "transform.h"
#include "debugging.h"
#include "psimrcc.h"
#include "updater.h"

namespace psi{ namespace psimrcc{

using namespace std;

/*!
 * Runs a MRPT2 and a MRCCSD computation
 * @todo move this code in the CCMRCC class
 */
void mrccsd(Options & options)
{
  // Initialize the mp2 module (integrals,fock matrix(ces),denominators)
  CCMRCC        mrcc(options);

  if(options.get_bool("PERT_CBS") && options.get_bool("PERT_CBS_COUPLING")){
    mrcc.compute_first_order_amps();
  }

  options.print();
  // Initialize the appropriate updater
  Updater* updater;
//  if(options_get_str("CORR_ANSATZ")=="SR")
//    updater = static_cast<Updater*>(new MkUpdater());
  if(options.get_str("CORR_ANSATZ")=="MK")
    updater = dynamic_cast<Updater*>(new MkUpdater(options));
  if(options.get_str("CORR_ANSATZ")=="BW")
    updater = dynamic_cast<Updater*>(new BWUpdater(options));

	// Compute the energy
  mrcc.compute_energy(updater);

  if(options.get_bool("PERT_CBS")){
    mrcc.perturbative_cbs();
  }

  delete updater;
}

/*!
 * Runs a CCSD_MP2 computation
 */
void mp2_ccsd(Options &options)
{
  // Initialize the mp2 module (integrals,fock matrix(ces),denominators)
  MP2_CCSD        mp2_ccsd(options);

  // Compute the initial amplitudes and CCSD_MP2 energy
  mp2_ccsd.compute_mp2_ccsd_energy();

  DEBUGGING(1,
    blas->print_memory();
  )
}

/*!
 * Runs a MRPT2 computation
 */
void mrpt2(Options &options)
{
  // Initialize the mp2 module (integrals,fock matrix(ces),denominators)
  IDMRPT2        idmrpt2(options);

  Updater* updater = dynamic_cast<Updater*>(new MkUpdater(options));

  // Compute the initial amplitudes and MP2 energy
  idmrpt2.compute_mrpt2_energy(updater);

  delete updater;

  DEBUGGING(1,
    blas->print_memory();
  )
}

/*!
 * Runs a integral transformation
 * @todo CCTransform is still unused in the code
 */
void transform_integrals()
{
//   CCTransform transf;
//   transf.read_so_integrals();
}

}} /* End Namespaces */
