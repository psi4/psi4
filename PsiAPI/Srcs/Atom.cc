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
#include<sstream>
#include<iomanip>
#include "Atom.h"
#include "libPsiUtil/AtomicData.h"
namespace PsiAPI{


Atom::Atom(size_t Z,const double* Carts, double Q, double Mass):
	Z_(Z),Carts_(Carts),Q_(Q){
	static const psi::AtomicData CRC;
	Symbol_ = (CRC[Z].AtSym());
	Mass_ = (Mass<0.0? CRC[Z].Mass() :Mass);
}

Atom::operator std::string()const{
   std::stringstream Result;
   Result<<Symbol_<<" "<<std::fixed
         <<std::setprecision(15)
         <<(*this)[0]<<" "<<(*this)[1]<<" "<<(*this)[2];
   return Result.str();
}


}


