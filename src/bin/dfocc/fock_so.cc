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
#include <fstream>
#include <psifiles.h>
#include <libiwl/iwl.hpp>
#include <libqt/qt.h>
#include "dfocc.h"

using namespace psi;
using namespace std;

namespace psi{ namespace dfoccwave{

void DFOCC::fock_so()
{      
  //psi::fprintf(outfile,"\tfock_so is starting... \n"); fflush(outfile);
/********************************************************************************************/
/************************** Build the Fock matrix *******************************************/
/********************************************************************************************/
     // Read From wfn
     Fa_ = SharedMatrix(reference_wavefunction_->Fa());
     //Fa_->print();

     // Memalloc
     FsoA = SharedTensor2d(new Tensor2d("SO-basis Alpha Fock Matrix", nso_, nso_));
     DsoA = SharedTensor2d(new Tensor2d("SO-basis Alpha Density Matrix", nso_, nso_));
     DQmatA = SharedTensor2d(new Tensor2d("Alpha D_munu^Q", nQ_ref, nso2_));
     DQvecA = SharedTensor1d(new Tensor1d("Alpha D^Q", nQ_ref));

     // Dmn = \sum_{i} Cmi Cni; Here, I used spin free density matrix instead of Alpha spin component
     DsoA->gemm(false, true, CoccA, CoccA, 1.0, 0.0);
     DsoA->scale(2.0);

     // D^Q = \sum_{mn} b_mn^Q Dmn
     DQvecA->gemv(false, nQ_ref, nso_ * nso_, bQso, DsoA, 1.0, 0.0);

     // D_mn^Q = \sum_{s} Dms * B_ns^Q
     DQmatA->contract233(false, false, nso_, nso_, DsoA, bQso, 1.0, 0.0);

     // 1e-part
     FsoA->copy(Hso);
 
     // Fmn += \sum_{Q} b_mn^Q D^Q
     FsoA->gemv(true, bQso, DQvecA, 1.0, 1.0);

     // Fmn += - 0.5 \sum_{Q,l} b_ml^Q D_nl^Q
     FsoA->contract(true, false, nso_, nso_, nso_ * nQ_ref, bQso, DQmatA, -0.5, 1.0);
     //FsoA->print();

//psi::fprintf(outfile,"\tfock_so is done. \n"); fflush(outfile);
}// end fock_so
}} // End Namespaces




