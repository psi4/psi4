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

#include "GMBE.h"
#include "MBEFrag.h"
#include "psi4-dec.h"

namespace psi{
namespace LibFrag {
double GMBE::Energy(const std::vector<MBEFragSet>& Systems,
      const std::vector<boost::shared_ptr<double[]> >& Energies,
      std::string& RealName) {
   double TEnergy=0.0;
   /*If N==1 special case and energies are only in Energies[0],but
   //we get sizes from Systems[1] and Systems[2]
   int index1=(N==1?0:1);
   for (int i=0; i<Systems[1].size(); i++)
      TEnergy+=Systemsts[i]*Energies[index1][i];

   for (int j=0; j<Systems[2].size(); j++) {
      int index2=(N==1?Systems[1].size()+j:j);
      TEnergy-=NMults[j]*Energies[index1][index2];
   }
   psi::outfile->Printf( "The total %d-body GMBE %s is: %16.12f (a.u.)", N,
         RealName.c_str(),TEnergy);*/
   return TEnergy;
}
}}//End namespaces

