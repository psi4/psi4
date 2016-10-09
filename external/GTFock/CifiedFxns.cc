/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2016 The Psi4 Developers.
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

#include <omp.h>
extern "C" {
#include "CifiedFxns.h"
}
#include "../libqt/qt.h"
#include "gtfock/libcint/basisset.h"
#include "../libmints/twobody.h"
#include "../libmints/integral.h"
#include "../libmints/basisset.h"

void timer_interface_on(char *name){
    psi::timer_on(name);
}
void timer_interface_off(char *name){
    psi::timer_off(name);
}

typedef psi::IntegralFactory IntFac;
typedef boost::shared_ptr<IntFac> SharedFac;
typedef boost::shared_ptr<psi::TwoBodyAOInt> SharedInts;
typedef boost::shared_ptr<psi::BasisSet> SharedBasis;
std::vector<SharedInts> Ints;
static SharedBasis primary_;
static SharedFac factory_;

void SetUp(){
    psi::Options& options = psi::Process::environment.options;
    primary_=psi::BasisSet::pyconstruct_orbital(
        psi::Process::environment.legacy_molecule(),
        "BASIS", options.get_str("BASIS")
    );
   factory_=SharedFac(new IntFac(primary_,primary_,primary_,primary_));
   // Each OMP thread gets its own integral object.
   // Make sure the vector is empty before pushing.
   Ints.clear();
   int maxthreads = omp_get_max_threads();
   for(int thread = 0; thread < maxthreads; ++thread) {
       Ints.push_back(SharedInts(factory_->erd_eri()));
   }
}

int ComputeShellQuartet(struct BasisSet*,int ThreadID,
      int M,int N,int P,int Q,double ** IntOut){
   int NInts=Ints[ThreadID]->compute_shell(M,N,P,Q);
   const double* Intergrals;
   if(NInts!=0)Intergrals=Ints[ThreadID]->buffer();
   (*IntOut)=const_cast<double *>(Intergrals);
   return NInts;
}

