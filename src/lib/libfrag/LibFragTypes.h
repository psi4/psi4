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
#ifndef LIBFRAGTYPES_H_
#define LIBFRAGTYPES_H_

#include <boost/shared_ptr.hpp>
#include <vector>
#include <sstream>
#include "physconst.h"
namespace psi{
   class Molecule;

namespace LibFrag{
///Typedefs I use all over the place
typedef psi::Molecule Mol;
typedef boost::shared_ptr<Mol> SharedMol;
class MBEFrag;
typedef boost::shared_ptr<MBEFrag> SharedFrag;
typedef std::vector<SharedFrag> NMerSet;

/*//Tiny wrapper class for atomic properties
class Atom{
   private:

      ///Function that actually does the copy
      void Copy(const Atom& other);

   public:

      ///The mass of the current atom
      double mass;

      ///The cartesian coordinates (a.u.) of the current atom
      std::vector<double> carts;

      ///Default constructor
      Atom():carts(3,0),mass(0){}

      ///Default destructor
      virtual ~Atom(){}

      ///Copy constructor
      Atom(const Atom& other){this->Copy(other);}

      ///Assignment operator
      const Atom& operator=(const Atom& other){
         if(this!=&other)this->Copy(other);return *this;
      }
};

class Cap:public Atom{
   private:
      std::string Symbol;
      int AtomIReplaced;
      void Copy(const Cap& other){
         this->Symbol=other.Symbol;
         this->AtomIReplaced=other.AtomIReplaced;
      }
   public:
      int ReplacedAtom(){return AtomIReplaced;}
      Cap(std::string Sym,double x,double y,double z,int ReplacedAtom):
         Atom(),AtomIReplaced(ReplacedAtom),Symbol(Sym){
         carts[0]=x;carts[1]=y;carts[2]=z;
      }
      Cap(const Cap& Other):Atom(Other){this->Copy(Other);}
      const Cap& operator=(const Cap& Other){
         if (this!=&Other) {
            Atom::operator=(Other);
            this->Copy(Other);
         }
         return *this;
      }
      std::string print_out(){
         std::stringstream temp;
         temp<<Symbol;
         for (int i=0; i<3; i++)temp<<" "<<carts[i]*pc_bohr2angstroms;
         temp<<"\n";
         return temp.str();
      }

};*/

}}//End namespaces



#endif /* LIBFRAGTYPES_H_ */
