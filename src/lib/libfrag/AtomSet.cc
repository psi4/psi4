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
#include "masses.h"

namespace psi {
namespace LibFrag {

/***********Atom Object Functions***********/
void Atom::Copy(const Atom& other) {
   this->Z_=other.Z_;
   this->mass_=other.mass_;
}

const Atom& Atom::operator=(const Atom& other) {
   CartObject::operator=(other);
   if (this!=&other) this->Copy(other);
   return *this;
}

Atom::Atom(const int Z, const double m, const double*Carts) :
      Z_(Z), mass_(m), CartObject(Carts) {
}

Atom::Atom(const int Z, const double m, const double x, const double y,
      const double z) :
      Z_(Z), mass_(m), CartObject(x, y, z) {
}

/*********Stand In Functions***************/

StandIn::StandIn(const int R, const double x, const double y, const double z) :
      CartObject(x, y, z), AtomIReplace_(R) {
}

StandIn::StandIn(const int R, const double* Carts) :
      CartObject(Carts), AtomIReplace_(R) {
}

const StandIn& StandIn::operator=(const StandIn& other) {
   CartObject::operator=(other);
   if (this!=&other) this->Copy(other);
   return *this;
}

/***********Cap Functions*******************/
Cap::Cap(const int BA, const int RA, const int Z, const int m,
      const double* Carts) :
      CartObject(Carts), StandIn(RA, Carts), Atom(Z, m, Carts),
            AtomIBonded2_(BA) {
}
Cap::Cap(const int BondAtom, const int ReplaceAtom, const int Z, const int m,
      const double x, const double y, const double z) :
      StandIn(ReplaceAtom, x, y, z), Atom(Z, m, x, y, z), CartObject(x, y, z),
            AtomIBonded2_(BondAtom) {
}

Cap::Cap(const Cap& other) :
      StandIn(other), Atom(other), CartObject(other) {
   this->Copy(other);
}

const Cap& Cap::operator=(const Cap& other) {
   Atom::operator=(other);
   StandIn::operator=(other);
   CartObject::operator=(other);
   if (this!=&other) this->Copy(other);
   return *this;
}

std::string Cap::print_out()const{
 std::stringstream message;
 message<<atomic_labels[this->Z()]<<" ";
 for(int i=0;i<3;i++)message<<this->Carts()[i]<<" ";
 message<<"\n";
 return message.str();
}
/**********Ghost Functions (Spooky!!)*********/

Ghost::Ghost(const int RA, const int Z, const double* Carts) :
      StandIn(RA, Carts), Atom(Z, 0.0, Carts), CartObject(Carts) {
}

Ghost::Ghost(const int RA, const int Z, const double x, const double y,
      const double z) :
      StandIn(RA, x, y, z), Atom(RA, 0.0, x, y, z), CartObject(x, y, z) {
}

Ghost::Ghost(const Ghost& other) :
      StandIn(other), Atom(other), CartObject(other) {
}

const Ghost& Ghost::operator=(const Ghost& other) {
   StandIn::operator=(other);
   Atom::operator=(other);
   CartObject::operator=(other);
   return *this;
}

}
} //End namespaces

