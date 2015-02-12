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
#ifndef SRC_LIB_LIBFRAG_LIBMOLECULE_PARAMETER_H_
#define SRC_LIB_LIBFRAG_LIBMOLECULE_PARAMETER_H_
#include <sstream>
namespace psi {
namespace LibMolecule {

class GeometricQuantity{
   public:
      ///Prints this quantity
      virtual std::string PrintOut()const=0;
      virtual ~GeometricQuantity(){}
      virtual int size()const=0;
      virtual int operator[](const int i)const=0;
};


template <int T>
class Parameter: public GeometricQuantity {
   private:
      int Members_[T];
      double Value_;
   protected:
      std::string Name_;
      std::string Unit_;
   public:
      int size()const{return T;}
      int operator[](const int i)const{return Members_[i];}
      Parameter(int Members[], double Value) :
            Value_(Value) {
         memcpy(Members_, Members, sizeof(int)*T);
      }
      double GetValue() const {
         return Value_;
      }
      std::string PrintOut() const {
         std::stringstream Message;
         Message<<Name_<<": ";
         for (int i=0; i<T; i++)
            Message<<(*this)[i]<<" ";
         Message<<"Value ("<<Unit_<<"):"<<Value_<<std::endl;
         return Message.str();
      }
};

class Pair:public Parameter<2> {
   public:
      Pair(int Members[], double Value) :
            Parameter(Members, Value) {
         Name_="Pair";
         Unit_="a.u.";
      }
};

class Bond:public Parameter<2> {
   public:
      Bond(int Members[], double Value) :
            Parameter(Members, Value) {
         Name_="Bond";
         Unit_="a.u.";
      }
};

class Angle:public Parameter<3> {
   public:
      Angle(int Members[], double Value) :
            Parameter(Members, Value) {
         Name_="Angle";
         Unit_="radians";
      }
};

class Torsion:public Parameter<4> {
   public:
      Torsion(int Members[], double Value) :
            Parameter(Members, Value) {
         Name_="Torsion";
         Unit_="radians";
      }
};

}} //End namespaces

#endif /* SRC_LIB_LIBFRAG_LIBMOLECULE_PARAMETER_H_ */
