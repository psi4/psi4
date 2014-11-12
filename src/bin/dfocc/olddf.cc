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

/** Standard library includes */
#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <sstream>
#include <fstream>
#include <string>
#include <iomanip>
#include <vector>

/** Required PSI3 includes */ 
#include <psifiles.h>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>
#include <libchkpt/chkpt.h>
#include <libpsio/psio.hpp>
#include <libchkpt/chkpt.hpp>
#include <libiwl/iwl.h>
#include <libqt/qt.h>
#include <libtrans/mospace.h>
#include <libtrans/integraltransform.h>

/** Required libmints includes */
#include <libmints/mints.h>
#include <libmints/factory.h>
#include <libmints/wavefunction.h>
#include<lib3index/dftensor.h>

#include "defines.h"
#include "dfocc.h"

using namespace boost;
using namespace psi;
using namespace std;


namespace psi{ namespace dfoccwave{
  
void DFOCC::df()
{   
  boost::shared_ptr<BasisSet> auxiliary = BasisSet::pyconstruct_auxiliary(reference_wavefunction_->molecule(),
            "DF_BASIS_CC", options_.get_str("DF_BASIS_CC"), "RIFIT", options_.get_str("BASIS"));

  boost::shared_ptr<DFTensor> DF (new DFTensor(reference_wavefunction_->basisset(), auxiliary, Ca_, noccA, nvirA, naoccA, navirA, options_));
  nQ = auxiliary->nbf(); // reads number of aux-basis functions
  //bQnn = boost::shared_ptr<Matrix>(new Matrix("B_munu^Q", nQ, nso2_));
  //bQnn = DF->Qso(); // reads b(Q|mu nu) where mu/nu is NOT packed 

  // Read MO basis intermediates
  // (pq|rs)  = \sum_{Q} (pq|Q) (Q|rs)
  //bQmo = DF->Qmo(); // reads (Q|pq) where p/q is NOT packed
  SharedMatrix tmpQso = DF->Qso(); // reads b(Q|ij)
  SharedMatrix tmpQoo = DF->Qoo(); // reads b(Q|ij)
  SharedMatrix tmpQov = DF->Qov(); // reads b(Q|ia)
  SharedMatrix tmpQvv = DF->Qvv(); // reads b(Q|ab)
  tmpQso->print();
  tmpQoo->print();
  tmpQov->print();
  tmpQvv->print();

} // end rhf
}} // End Namespaces



