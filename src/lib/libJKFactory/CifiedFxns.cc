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
SharedInts Ints;

int ComputeShellQuartet(struct BasisSet*,int ThreadID,
      int M,int N,int P,int Q,double ** IntOut){
   psi::Options& options = psi::Process::environment.options;
   SharedBasis primary =psi::BasisSet::pyconstruct_orbital(
         psi::Process::environment.molecule(),
         "BASIS", options.get_str("BASIS")
   );
   SharedFac factory(new IntFac(primary,primary,primary,primary));
   Ints=SharedInts(factory->erd_eri());
   std::cout<<"Computing shell: "<<M<<" "<<N<<" "<<P<<" "<<Q<<std::endl;
   int NInts=Ints->compute_shell(M,N,P,Q);
   const double* Intergrals;
   if(NInts!=0)Intergrals=Ints->buffer();
   (*IntOut)=const_cast<double *>(Intergrals);
   return NInts;
}

