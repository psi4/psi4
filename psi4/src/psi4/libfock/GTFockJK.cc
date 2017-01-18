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
#include "psi4/libpsi4util/exception.h"
#include "psi4/libfock/jk.h"
#ifdef ENABLE_GTFOCK
#include <GTFock/MinimalInterface.h>
#else
namespace psi {
struct MinimalInterface{
    MinimalInterface(size_t,bool)
    {
        throw PSIEXCEPTION("PSI4 has not been compiled with GTFock support");
    }
    void SetP(std::vector<SharedMatrix>&){}
    void GetJ(std::vector<SharedMatrix>&){}
    void GetK(std::vector<SharedMatrix>&){}

};
}
#endif


#ifdef ENABLE_GTFOCK
namespace psi {
GTFockJK::GTFockJK(std::shared_ptr<psi::BasisSet> Primary) :
      JK(Primary),Impl_(new MinimalInterface()){

}
void GTFockJK::compute_JK() {

   NMats_ = C_left_.size();
   Impl_->create_pfock(NMats_, lr_symmetric_);
   Impl_->SetP(D_ao_);
   Impl_->GetJ(J_ao_);
   Impl_->GetK(K_ao_);
   Impl_->destroy_gtfock();
}


}
#endif
