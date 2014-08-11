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

#include "AtomSet.h"
#include <cmath>
namespace LibFrag{

void Atom::Copy(const Atom& other){
   this->mass=other.mass;
   this->carts=other.carts;
}


const AtomSet& AtomSet::operator=(const AtomSet& other) {
   if (this!=&other) {
      this->Set::operator=(other);
      this->Copy(other);
   }
   return *this;
}

void AtomSet::AddMass(const int i, const double m){
   CoM.clear();
   Elem2Atoms[i].mass=m;
}

void AtomSet::AddCarts(const int i, const double x, const double y,
      const double z) {
   CoM.clear();
   Elem2Atoms[i].carts[0]=x;
   Elem2Atoms[i].carts[1]=y;
   Elem2Atoms[i].carts[2]=z;
}

void AtomSet::CalcCoM() {
   double M=0.0;
   CoM.clear();
   for(int i=0;i<3;i++)CoM.push_back(0.0);
   for (int i=0; i<Atoms.size(); i++) {
      M+=Mass(i);
      std::vector<double> tempcarts=Carts(i);
      for (int j=0; j<3; j++)
         CoM[j]+=Mass(i)*tempcarts[j];
   }
   for (int i=0; i<3; i++)
      CoM[i]/=M;
}

double AtomSet::Distance(AtomSet&other) {
   if (this->CoM.size()!=3) this->CalcCoM();
   if (other.CoM.size()!=3) other.CalcCoM();
   double distsq=0.0;
   for (int i=0; i<3; i++)
      distsq+=(this->CoM[i]-other.CoM[i])*(this->CoM[i]-other.CoM[i]);
   return sqrt(distsq);
}

void AtomSet::Copy(const AtomSet& other) {
   this->Ghosts=other.Ghosts;
   this->Charges=other.Charges;
   this->Caps=other.Caps;
   this->Elem2Atoms=other.Elem2Atoms;
   this->CoM=other.CoM;
}

}//End namespace

