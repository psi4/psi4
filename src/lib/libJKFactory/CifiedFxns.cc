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
namespace Interface{
std::vector<boost::shared_ptr<psi::TwoBodyAOInt> > Ints(omp_get_max_threads());
}
bool IntsSet=false;

void SetUpInts(){
psi::Options& options = psi::Process::environment.options;
boost::shared_ptr<psi::BasisSet> primary =
		psi::BasisSet::pyconstruct_orbital(psi::Process::environment.molecule(),
    "BASIS", options.get_str("BASIS"));
boost::shared_ptr<psi::IntegralFactory>
factory(new psi::IntegralFactory(primary,primary,primary,primary));
for(int i=0;i<omp_get_max_threads();i++){
  Interface::Ints[i]=
     boost::shared_ptr<psi::TwoBodyAOInt>(factory->erd_eri());
}
}


void timer_interface_on(char *name){
    psi::timer_on(name);
}
void timer_interface_off(char *name){
    psi::timer_off(name);
}

int ComputeShellQuartet(struct BasisSet*,int NThreads,int M,int N,int P,int Q,
         double **IntsOut){
   if(!IntsSet){
	   SetUpInts();
	   IntsSet=true;
   }
	int NInts=Interface::Ints[NThreads]->compute_shell(M,N,P,Q);
   *IntsOut=
      const_cast<double *>(Interface::Ints[NThreads]->buffer());
   return NInts;
}

