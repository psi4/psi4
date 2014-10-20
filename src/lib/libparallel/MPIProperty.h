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
#ifndef MPIPROPERTY_H_
#define MPIPROPERTY_H_

namespace psi{

template<typename T>
class MPIProperty;

template<typename T>
class MPIProperty<T*>{
   protected:
      T* Property_;
      std::string Name_;
      int Length_;
   public:
      MPIProperty<T>(const T* Property,const int Length):
      Property_(Property),Length_(Length){}
      T GetProp()const{return *Property_;}
      int GetLength()const{return Length_;}
};

}



#endif /* MPIPROPERTY_H_ */
