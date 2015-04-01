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
#ifndef SRC_LIB_LIBFRAG_LIBMOLECULE_SUPERCELL_H_
#define SRC_LIB_LIBFRAG_LIBMOLECULE_SUPERCELL_H_

#include <boost/shared_ptr.hpp>
#include "UnitCell.h"

namespace psi{
namespace LibMolecule{
/** \brief A class for managing a unitcell replicated in each direction
 *
 *  Of particular usefulness this class keeps track of which atoms
 *  were in the original cell.  It's also worth knowing that
 *  the first cell of atoms contained within the supercell is the unitcell.
 *  I.e. if the unitcell contains 20 atoms, the first 20 atoms of the
 *  supercell returned by a MolItr are the 20 in the unitcell.
 */
class SuperCell: public UnitCell{
   private:
      boost::shared_ptr<UnitCell> OrigCell_;
      void FormSuperCell(const int dims[],const double cellvecs[],
            const UnitCell& UC, int depth=0);
   public:
      ///Replicates the Original Cell, NShells times in each direction
      SuperCell(const UnitCell& OrigCell, const std::vector<int>& NShells);
      boost::shared_ptr<const UnitCell> GetUnitCell()const{return OrigCell_;}
};


}}//End namespaces



#endif /* SRC_LIB_LIBFRAG_LIBMOLECULE_SUPERCELL_H_ */
