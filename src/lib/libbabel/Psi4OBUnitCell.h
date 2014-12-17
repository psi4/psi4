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
#ifndef PSI4OBUNITCELL_H_
#define PSI4OBUNITCELL_H_
#include "openbabel/mol.h"
#include "openbabel/generic.h"

namespace psi{

/** \brief Class created by Professor Rafal Podeszwa to avoid dangling bonds
 *               that occur in the normal OpenBabel OBUnitCell
 *
 *
 *   This class implements a FillUnitCell2 function that should be called
 *   instead of OpenBabel::OBUnitCell::FillUnitCell because
 *
 */
class OBUnitCellChild : public OpenBabel::OBUnitCell{
 public:
   void FillUnitCell2(OpenBabel::OBMol *mol);
   OBUnitCellChild(const OBUnitCell& other):OBUnitCell(other){}
};

}



#endif /* PSI4OBUNITCELL_H_ */
