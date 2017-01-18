/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
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
 * @END LICENSE
 */

#include"ccsd.h"
#include"frozen_natural_orbitals.h"
#include "psi4/libciomr/libciomr.h"
#include "psi4/libtrans/integraltransform.h"
#include "psi4/libtrans/mospace.h"

namespace psi{ namespace fnocc{

SharedWavefunction fnocc(SharedWavefunction ref_wfn, Options &options) {

  std::shared_ptr<Wavefunction> wfn;

  if ( !options.get_bool("DFCC") ){

      // frozen natural orbital ccsd(t)
      if (options.get_bool("NAT_ORBS")) {

          std::shared_ptr<FrozenNO> fno(new FrozenNO(ref_wfn, options));
          fno->ComputeNaturalOrbitals();
          wfn = (std::shared_ptr<Wavefunction>)fno;

      }else {
          wfn = ref_wfn;
      }

      // transform integrals
      tstart();
      outfile->Printf("        ==> Transform all two-electron integrals <==\n");
      outfile->Printf("\n");

      std::vector<std::shared_ptr<MOSpace> > spaces;
      spaces.push_back(MOSpace::all);
      std::shared_ptr<IntegralTransform> ints(new IntegralTransform(wfn, spaces, IntegralTransform::Restricted,
                 IntegralTransform::IWLOnly, IntegralTransform::QTOrder, IntegralTransform::OccAndVir, false));
      ints->set_dpd_id(0);
      ints->set_keep_iwl_so_ints(true);
      ints->set_keep_dpd_so_ints(true);
      ints->initialize();
      ints->transform_tei(MOSpace::all, MOSpace::all, MOSpace::all, MOSpace::all);
      tstop();

      if ( !options.get_bool("RUN_CEPA") ) {
          std::shared_ptr<CoupledCluster> ccsd(new CoupledCluster(wfn,options));
          ccsd->compute_energy();
      } else {
          std::shared_ptr<CoupledPair> cepa (new CoupledPair(wfn,options));
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
      std::shared_ptr<DFFrozenNO> fno(new DFFrozenNO(ref_wfn,options));
      fno->ThreeIndexIntegrals();
      if ( options.get_bool("NAT_ORBS") ) {
          fno->ComputeNaturalOrbitals();
          wfn = (std::shared_ptr<Wavefunction>)fno;
      }else {
          wfn = ref_wfn;
      }
      // ccsd(t)!

      #ifdef GPUCC
          std::shared_ptr<GPUDFCoupledCluster> ccsd (new GPUDFCoupledCluster(wfn,options));
          ccsd->compute_energy();
      #else
          std::shared_ptr<DFCoupledCluster> ccsd (new DFCoupledCluster(wfn,options));
          ccsd->compute_energy();
      #endif

      tstop();


  }

  return wfn;
} // end fnocc

}} // end namespaces
