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

#include "MBEFrag.h"

namespace psi {
namespace LibFrag {

void MBEFrag::Copy(const MBEFrag& other) {
   this->Parents=other.Parents;
   this->MBEOrder=other.MBEOrder;
   this->Atoms_=other.Atoms_;
   this->Caps_=other.Caps_;
   this->Charges_=other.Charges_;
   this->Ghosts_=other.Ghosts_;
}

void UpdateCaps(CartSet<SharedCap>& ThisCaps,
      const CartSet<SharedCap>& OtherCaps, MBEFrag* Frag) {
   ThisCaps*=OtherCaps;
   if (ThisCaps.size()>0) {
      CartSet<SharedCap> temp(ThisCaps);
      temp.Clear();
      //less any that are now in Atoms_
      const CartSet<SharedAtom> Atoms=Frag->Atoms();
      for (int i=0; i<ThisCaps.size(); i++) {
         int Atom2Check=ThisCaps.Object(i)->ReplacedAtom();
         if (!Atoms.Contains(Atom2Check)) {
            temp<<Atom2Check;
         }
      }
      ThisCaps=temp;
   }
}

void MBEFrag::operator*=(const MBEFrag& other) {
   this->Atoms_*=other.Atoms_;

   //For the Charges and and Ghosts need to subtract out new Atoms
   this->Charges_-=other.Atoms_;
   this->Ghosts_-=other.Atoms_;

   //For the Caps we have their union
   UpdateCaps(this->Caps_, other.Caps_, this);

}

void MBEFrag::operator/=(const MBEFrag& other) {
   //For the Charges and Ghosts need union with old Atoms_
   this->Ghosts_*=this->Atoms_;
   this->Charges_*=this->Atoms_;

   //Now take the intersection, to get our new set of atoms
   this->Atoms_/=other.Atoms_;

   //Subtract it out of the Charges and Ghosts
   this->Ghosts_-=this->Atoms_;
   this->Charges_-=this->Atoms_;
   UpdateCaps(this->Caps_, other.Caps_, this);
}

}
} //End namespaces

