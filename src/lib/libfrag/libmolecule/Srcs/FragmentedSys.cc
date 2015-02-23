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
#include <sstream>
#include "FragmentedSys.h"
#include "Utils/Fragmenter.h"
#include "LibFragFragment.h"
#include "PsiMap.h"
#include "MoleculeTypes.h"
#include "Utils/GeomManipulator.h"
namespace psi{
namespace LibMolecule{

typedef boost::shared_ptr<Fragment> SharedFrag;
typedef std::vector<SharedFrag> vSharedFrag;

int FragmentedSystem::N()const{return NMers_.N();}
double FragmentedSystem::Coef(const int N, const int i)const{
   return NMers_.ScaleFacts_(N,i);
}
FragmentedSystem::FragmentedSystem(const Molecule& System2Frag,
      const int N):FragSysGuts(System2Frag,N){}
FragmentedSystem::FragmentedSystem(const SuperCell& System2Frag,
      const int N):FragSysGuts(System2Frag,N){}

std::string FragmentedSystem::PrintOut(const int Value)const{
   std::stringstream Result;
   Result<<NMers_.PrintOut(Value);
   Result<<NMers_.ScaleFacts_.PrintOut()<<std::endl;
   return Result.str();
}

FragmentedSystem::iterator FragmentedSystem::begin(const int N)const{
   return NMers_.NMerBegin(N+1);
}

FragmentedSystem::iterator FragmentedSystem::end(const int N)const{
   return NMers_.NMerEnd(N+1);
}

}}//End namespaces

