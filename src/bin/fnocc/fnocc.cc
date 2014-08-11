/*
 *@BEGIN LICENSE
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

#include"ccsd.h"
#include"frozen_natural_orbitals.h"
#include<libciomr/libciomr.h>
#include<libtrans/integraltransform.h>
#include<libtrans/mospace.h>

using namespace boost;

namespace psi{ namespace fnocc{

PsiReturnType fnocc(Options &options) {

  boost::shared_ptr<Wavefunction> wfn;

  if ( !options.get_bool("DFCC") ){

      // frozen natural orbital ccsd(t)
      if (options.get_bool("NAT_ORBS")) {

          boost::shared_ptr<FrozenNO> fno(new FrozenNO(Process::environment.wavefunction(),options));
          fno->ComputeNaturalOrbitals();
          wfn = (boost::shared_ptr<Wavefunction>)fno;

          // transform integrals
          tstart();
          outfile->Printf("        ==> Transform all two-electron integrals <==\n");
          outfile->Printf("\n");

          std::vector<shared_ptr<MOSpace> > spaces;
          spaces.push_back(MOSpace::all);
          boost::shared_ptr<IntegralTransform> ints(new IntegralTransform(wfn, spaces, IntegralTransform::Restricted,
                     IntegralTransform::IWLOnly, IntegralTransform::QTOrder, IntegralTransform::OccAndVir, false));
          ints->set_dpd_id(0);
          ints->set_keep_iwl_so_ints(true);
          ints->set_keep_dpd_so_ints(true);
          ints->initialize();
          ints->transform_tei(MOSpace::all, MOSpace::all, MOSpace::all, MOSpace::all);
          tstop();

      }else {
          wfn = Process::environment.wavefunction();
      }

      if ( !options.get_bool("RUN_CEPA") ) {
          boost::shared_ptr<CoupledCluster> ccsd(new CoupledCluster(wfn,options));
          Process::environment.set_wavefunction(ccsd);
          ccsd->compute_energy();
      } else {
          boost::shared_ptr<CoupledPair> cepa (new CoupledPair(wfn,options));
          Process::environment.set_wavefunction(cepa);
          cepa->compute_energy();
      }

  }else {

      tstart();
      
      outfile->Printf("\n\n");
      outfile->Printf( "        *******************************************************\n");
      outfile->Printf( "        *                                                     *\n");
      outfile->Printf( "        *                       DF-CCSD                       *\n");
      outfile->Printf( "        *                 Density-fitted CCSD                 *\n");
      outfile->Printf( "        *                                                     *\n");
      outfile->Printf( "        *                   Eugene DePrince                   *\n");
      outfile->Printf( "        *                                                     *\n");
      outfile->Printf( "        *******************************************************\n");
      outfile->Printf("\n\n");
      

      // three-index integrals are generated/read by fno class
      boost::shared_ptr<DFFrozenNO> fno(new DFFrozenNO(Process::environment.wavefunction(),options));
      fno->ThreeIndexIntegrals();
      if ( options.get_bool("NAT_ORBS") ) {
          fno->ComputeNaturalOrbitals();
          wfn = (boost::shared_ptr<Wavefunction>)fno;
      }else {
          wfn = Process::environment.wavefunction();
      }
      // ccsd(t)!
      boost::shared_ptr<DFCoupledCluster> ccsd (new DFCoupledCluster(wfn,options));
      ccsd->compute_energy();
      tstop();
  }

  return  Success;
} // end fnocc

}} // end namespaces
