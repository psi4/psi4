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
#include "Geometry.h"
#include "LibFragMolecule.h"
#include "AtomicData.h"
namespace psi {
namespace LibMolecule {
typedef boost::shared_ptr<Bond> SharedBond;
typedef boost::shared_ptr<Pair> SharedPair;
void Geometry::FormDistance() {
   MolItr Start=Mol_.Begin(),Done=Mol_.End();
   int pair[2];
   AtomicData CRC;
   for (int i=0; Start!=Done; ++Start, ++i) {
      double radius1=CRC[(*Start)->Z()].CovRad();
      pair[0]=i;
      double dr[3];
      for (int j=0; j<3; j++)
         dr[j]=(*(*Start))[j];
      MolItr Atom2=Start;
      for (int j=i; Atom2!=Done; ++Atom2, j++) {
         if (i==j) continue;
         double radius2=CRC[(*Atom2)->Z()].CovRad();
         pair[1]=j;
         double d2=0;
         for (int k=0; k<3; k++) {
            double dist=(dr[k]-(*(*Atom2))[k]);
            d2+=dist*dist;
         }
         double d=sqrt(d2);
         if (d<=1.2*(radius1+radius2)) {
            Bonds_.push_back(SharedBond(new Bond(pair, d)));
            Connections_[pair[0]].push_back(pair[1]);
            Connections_[pair[1]].push_back(pair[0]);
         }
         else Pairs_.push_back(SharedPair(new Pair(pair, d)));
      }
   }
}

std::string Connections::PrintOut()const{
   std::stringstream Result;
   for(int i=0;i<(*this).size();i++){
      for(int j=0;j<(*this)[i].size();j++)
         Result<<(*this)[i][j]<<" ";
      Result<<std::endl;
   }
   return Result.str();
}

Geometry::Geometry(const Molecule& Mol) :
      Distance_(new double[Mol.NAtoms()*Mol.NAtoms()]),
            Connections_(Mol.NAtoms()), Mol_(Mol){
   FormDistance();
}

}
} //End namespaces
