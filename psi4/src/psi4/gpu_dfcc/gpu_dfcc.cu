/*
 *@BEGIN LICENSE
 *
 * GPU-accelerated density-fitted coupled-cluster, a plugin to:
 *
 * PSI4: an ab initio quantum chemistry software package
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 *@END LICENSE
 */

#include <psi4/libplugin/plugin.h>
#include <psi4/psi4-dec.h>
#include "psi4/liboptions/liboptions.h"
#include <psi4/libpsio/psio.hpp>
#include "dfcc_fnocc.h"
#include "frozen_natural_orbitals.h"
#include "ccsd.h"


INIT_PLUGIN

using namespace std;
namespace psi{ namespace fnocc {

#if 0
extern "C" 
int read_options(std::string name, Options& options)
{
    if (name == "GPU_DFCC"|| options.read_globals()) {
      /*- Do dgemm timings? ! expert -*/
      options.add_bool("DGEMM_TIMINGS",false);
      /*- Maximum amount of pinned CPU memory (mb) -*/
      options.add_int("MAX_MAPPED_MEMORY",7000);
      /*- Override number of GPUs detected? -*/
      options.add_int("NUM_GPUS",0);
      /*- Do time each cc diagram? -*/
      options.add_bool("CC_TIMINGS",false);
      /*- Convergence for the CC energy.  Note that convergence is
          met only when E_CONVERGENCE and R_CONVERGENCE are satisfied. -*/
      options.add_double("E_CONVERGENCE", 1.0e-8);
      /*- Convergence for the CC amplitudes.  Note that convergence is
          met only when E_CONVERGENCE and R_CONVERGENCE are satisfied. -*/
      options.add_double("R_CONVERGENCE", 1.0e-7);
      /*- Maximum number of CC iterations -*/
      options.add_int("MAXITER", 100);
      /*- Desired number of DIIS vectors -*/
      options.add_int("DIIS_MAX_VECS", 8);
      /*- Do use low memory option for triples contribution? Note that this 
          option is enabled automatically if the memory requirements of the 
          conventional algorithm would exceed the available resources -*/
      //options.add_bool("TRIPLES_LOW_MEMORY",false);
      /*- Do compute triples contribution? !expert -*/
      options.add_bool("COMPUTE_TRIPLES", true);
      /*- Do use MP2 NOs to truncate virtual space for CCSD and (T)? -*/
      options.add_bool("NAT_ORBS", false);
      /*- Cutoff for occupation of MP2 NO orbitals in FNO-CCSD(T)
          ( only valid if |gpu_dfcc__nat_orbs| = true ) -*/
      options.add_double("OCC_TOLERANCE", 1.0e-6);
      /*- Do SCS-MP2? -*/
      options.add_bool("SCS_MP2", false);
      /*- Do SCS-CCSD? -*/
      options.add_bool("SCS_CCSD", false);
      /*- Do SCS-CEPA? Note that the scaling factors will be identical
      to those for SCS-CCSD. -*/
      options.add_bool("SCS_CEPA", false);
      /*- Opposite-spin scaling factor for SCS-MP2 -*/
      options.add_double("MP2_SCALE_OS",1.20);
      /*- Same-spin scaling factor for SCS-MP2 -*/
      options.add_double("MP2_SCALE_SS",1.0/3.0);
      /*- Oppposite-spin scaling factor for SCS-CCSD -*/
      options.add_double("CC_SCALE_OS", 1.27);
      /*- Same-spin scaling factor for SCS-CCSD -*/
      options.add_double("CC_SCALE_SS",1.13);

      /*- Do use density fitting in CC? This keyword is used internally
          by the driver. Changing its value will have no effect on the 
          computation. ! expert -*/
      options.add_bool("DFCC",false);
      /*- Auxilliary basis for df-ccsd(t). -*/
      options.add_str("DF_BASIS_CC","");
      /*- Tolerance for Cholesky decomposition of the ERI tensor. 
          ( only valid if |gpu_dfcc__df_basis_cc| = cholesky or 
          |scf__scf_type|=cd -*/
      options.add_double("CHOLESKY_TOLERANCE",1.0e-4);

    }

    return true;
}
#endif

SharedWavefunction gpu_dfcc(SharedWavefunction ref_wfn, Options &options) {

    std::shared_ptr<Wavefunction> wfn;
    std::shared_ptr<DFFrozenNO> fno(new DFFrozenNO(ref_wfn,options));
    fno->ThreeIndexIntegrals();

    if ( options.get_bool("NAT_ORBS") ) {
        fno->ComputeNaturalOrbitals();
        wfn = (std::shared_ptr<Wavefunction>)fno;
    }else {
        wfn = ref_wfn;
    }
    std::shared_ptr<GPUDFCoupledCluster> ccsd (new GPUDFCoupledCluster(wfn,options));
    ccsd->compute_energy();
    //return wfn;
    return (SharedWavefunction)ccsd;
}

}} // End namespaces

