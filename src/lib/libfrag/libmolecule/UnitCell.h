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
#ifndef SRC_LIB_LIBFRAG_LIBMOLECULE_UNITCELL_H_
#define SRC_LIB_LIBFRAG_LIBMOLECULE_UNITCELL_H_
#include <vector>
#include "Implementations/UnitCellGuts.h"
namespace psi{
namespace LibMolecule{
class Molecule;
class UnitCell: public UnitCellGuts{
   protected:
      ///This is only to be called by SuperCell, which will do it's own setup
      UnitCell(){}
   public:
      ///Returns an array containing the three angles (in Radians)
      const double* Angles()const{return UnitCellGuts::Angles();}
      ///Returns an array containing the three sides (in a.u.)
      const double* Sides()const{return UnitCellGuts::Sides();}
      /** \brief Creates a unit cell with the atoms in Mol
       *
       *  \param[in] The molecule that forms the unit cell
       *  \param[in] Sides The lattice vectors in the order a,b,c
       *  \param[in] Angles The lattice angles in the order alpha,beta,gamma
       *  \param[in] IsFrac True if the coordinates contained in the molecule
       *             are fractional coordinates
       *  \param[in] IsBhor True if the sides given are in atomic units
       *  \param[in] IsDegree True if the angles given are in degrees
       */
      UnitCell(const Molecule& Mol,const double* Sides,
            const double* Angles, const bool IsFrac=true,
            const bool IsBohr=false,const bool IsDegree=true):
               UnitCellGuts(Mol,Sides,Angles,IsFrac,IsBohr,IsDegree){}
     UnitCell(const UnitCell& other):UnitCellGuts(other){}
     void FixUnitCell();
     const UnitCell& operator=(const UnitCell& other){
        UnitCellGuts::operator=(other);
        return *this;
     }
};

}}



#endif /* SRC_LIB_LIBFRAG_LIBMOLECULE_UNITCELL_H_ */
